#!/usr/bin/python3
import sys
sys.path.append('../IrisUtils/')

import SoapySDR
from SoapySDR import * #SOAPY_SDR_ constants
from optparse import OptionParser
import numpy as np
import time
import os
import signal
import math
import pdb 
import json
import scipy.io as sio 
from functools import partial
from scipy.linalg import hadamard
from type_conv import *
from file_rdwr import *
from generate_sequence import *

sdrs = None
hub_dev = None

def beamsweeper(hub, serials, rate, freq, txgain, rxgain, numSamps, numSyms, prefix_length, postfix_length, calibrate, both_channels):
    global sdrs, hub_dev 
    if hub != "": hub_dev = SoapySDR.Device(dict(serial=hub))
    print("setting %s as eNB" % (serials))
    sdrs = [SoapySDR.Device(dict(serial=serial1)) for serial1 in serials.split(',')]

    #some default sample rates
    for sdr in sdrs:
        info = sdr.getHardwareInfo();
        print("%s settings on device" % (info["frontend"]))
        for ch in [0,1]:
            sdr.setBandwidth(SOAPY_SDR_RX, ch, 2.5*rate)
            sdr.setBandwidth(SOAPY_SDR_TX, ch, 2.5*rate)
            sdr.setSampleRate(SOAPY_SDR_TX, ch, rate)
            sdr.setSampleRate(SOAPY_SDR_RX, ch, rate)
            sdr.setFrequency(SOAPY_SDR_TX, ch, "BB", 0.75*rate)
            sdr.setFrequency(SOAPY_SDR_RX, ch, "BB", 0.75*rate)
            sdr.setFrequency(SOAPY_SDR_TX, ch, "RF", freq - 0.75*rate)
            sdr.setFrequency(SOAPY_SDR_RX, ch, "RF", freq - 0.75*rate)
            if ("CBRS" in info["frontend"]):
                sdr.setGain(SOAPY_SDR_TX, ch, 'ATTN', 0) #[-18,0] by 3
                sdr.setGain(SOAPY_SDR_TX, ch, 'PA1', 15) #[0,15]
                sdr.setGain(SOAPY_SDR_TX, ch, 'PA2', 0) #[0,15]
                sdr.setGain(SOAPY_SDR_TX, ch, 'PA3', 30) #[0,30]
            sdr.setGain(SOAPY_SDR_TX, ch, 'IAMP', 0) #[0,12]
            sdr.setGain(SOAPY_SDR_TX, ch, 'PAD', txgain) #[-52,0]

            if ("CBRS" in info["frontend"]):
                sdr.setGain(SOAPY_SDR_RX, ch, 'ATTN', 0) #[-18,0]
                sdr.setGain(SOAPY_SDR_RX, ch, 'LNA1', 30) #[0,33]
                sdr.setGain(SOAPY_SDR_RX, ch, 'LNA2', 17) #[0,17]
            sdr.setGain(SOAPY_SDR_RX, ch, 'LNA', rxgain) #[0,30]
            sdr.setGain(SOAPY_SDR_RX, ch, 'TIA', 0) #[0,12]
            sdr.setGain(SOAPY_SDR_RX, ch, 'PGA', 0) #[-12,19]
            sdr.setAntenna(SOAPY_SDR_RX, ch, "TRX")
        for ch in [0,1]:
            if calibrate:
                sdr.writeSetting(SOAPY_SDR_RX, ch, "CALIBRATE", 'SKLK')
                sdr.writeSetting(SOAPY_SDR_TX, ch, "CALIBRATE", '')
            sdr.setDCOffsetMode(SOAPY_SDR_RX, ch, True)

        if not both_channels:
            sdr.writeSetting(SOAPY_SDR_RX, 1, 'ENABLE_CHANNEL', 'false')
            sdr.writeSetting(SOAPY_SDR_TX, 1, 'ENABLE_CHANNEL', 'false')
        sdr.writeSetting("RESET_DATA_LOGIC", "")

    if hub == "":
        sdrs[0].writeSetting("SYNC_DELAYS", "")
    else:
        hub_dev.writeSetting("SYNC_DELAYS", "")

    #packet size
    symSamp = numSamps + prefix_length + postfix_length
    print("numSamps = %d" % symSamp)

    # preambles to be sent from BS and correlated against in UE
    upsample = 1
    preambles_bs = generate_training_seq(preamble_type='gold_ifft', seq_length=128, cp=0, upsample=1)
    preambles = preambles_bs[:,::upsample] #the correlators can run at lower rates, so we only need the downsampled signal.
    beacon = preambles[0,:]*.25

    possible_dim = []
    nRadios = len(sdrs)
    nChannels = 2 if both_channels else 1
    numAnt = nRadios * nChannels
    possible_dim.append(2**(np.ceil(np.log2(numAnt))))
    h_dim = min(possible_dim)
    hadamard_matrix = hadamard(h_dim)       #hadamard matrix : http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.linalg.hadamard.html
    beacon_weights = hadamard_matrix[0:numAnt, 0:numAnt]
    print(beacon_weights)
    beacon_weights = beacon_weights.astype(np.uint32)
    bzeros = np.array([0]*numAnt, np.uint32) 

    bsched = "BG"#+''.join("G"*(numSyms-1)) 
    print("Schedule %s " % bsched) 
    bconf = {"tdd_enabled": True,
             "frame_mode":
             "free_running",
             "symbol_size" : symSamp,
             "frames": [bsched],
             "beacon_start" : prefix_length,
             "beacon_stop" : prefix_length+len(beacon)}
    for i, sdr in enumerate(sdrs):
        sdr.writeSetting("TDD_CONFIG", json.dumps(bconf))
        sdr.writeSetting("TX_SW_DELAY", str(30))
        sdr.writeSetting("TDD_MODE", "true")

    for i, sdr in enumerate(sdrs):
        sdr.writeRegisters("BEACON_RAM", 0, cfloat2uint32(beacon, order='QI').tolist())
        sdr.writeRegisters("BEACON_RAM_WGT_A", 0, beacon_weights[i*nChannels].tolist())
        sdr.writeRegisters("BEACON_RAM_WGT_B", 0, beacon_weights[2*i+1].tolist() if both_channels else bzeros.tolist())
        sdr.writeSetting("BEACON_START", str(numAnt))

    signal.signal(signal.SIGINT, partial(signal_handler, rate, numSyms))
    if hub == "":
        sdrs[0].writeSetting("TRIGGER_GEN", "")
    else:
        hub_dev.writeSetting("TRIGGER_GEN", "")
    signal.pause()

def signal_handler(rate, numSyms, signal, frame):
    global sdrs, hub_dev
    print("printing number of frames")
    for sdr in sdrs:
        print("NB 0x%X" % SoapySDR.timeNsToTicks(sdr.getHardwareTime(""), rate))
        # ADC_rst, stops the tdd time counters
        sdr.writeSetting("RESET_DATA_LOGIC", "")
        conf = {"tdd_enabled" : False}
        sdr.writeSetting("TDD_CONFIG", json.dumps(conf))
        sdr.writeSetting("TDD_MODE", "false")
    sdrs = None
    hub_dev = None
    sys.exit(0)

def main():
    parser = OptionParser()
    parser.add_option("--hub", type="string", dest="hub", help="serial number of the hub device", default="")
    parser.add_option("--serials", type="string", dest="serials", help="serial numbers of the devices", default="")
    parser.add_option("--rate", type="float", dest="rate", help="Tx sample rate", default=5e6)
    parser.add_option("--txgain", type="float", dest="txgain", help="Optional Tx gain (dB)", default=40.0)
    parser.add_option("--rxgain", type="float", dest="rxgain", help="Optional Rx gain (dB)", default=20.0)
    parser.add_option("--freq", type="float", dest="freq", help="Optional Tx freq (Hz)", default=3.6e9)
    parser.add_option("--numSamps", type="int", dest="numSamps", help="Num samples to receive", default=512)
    parser.add_option("--prefix-length", type="int", dest="prefix_length", help="prefix padding length for beacon and pilot", default=82)
    parser.add_option("--postfix-length", type="int", dest="postfix_length", help="postfix padding length for beacon and pilot", default=68)
    parser.add_option("--numSyms", type="int", dest="numSyms", help="Number of symbols in one frame", default=20)
    parser.add_option("--calibrate", action="store_true", dest="calibrate", help="transmit from both channels",default=False)
    parser.add_option("--both-channels", action="store_true", dest="both_channels", help="transmit from both channels",default=False)
    (options, args) = parser.parse_args()
    beamsweeper(
        hub=options.hub,
        serials=options.serials,
        rate=options.rate,
        freq=options.freq,
        txgain=options.txgain,
        rxgain=options.rxgain,
        numSamps=options.numSamps,
        numSyms=options.numSyms,
        prefix_length=options.prefix_length,
        postfix_length=options.postfix_length,
        calibrate=options.calibrate,
        both_channels=options.both_channels
    )

if __name__ == '__main__':
    main()
 