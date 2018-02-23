#!/usr/bin/env python

import sys
import json
import ROOT
from optparse import OptionParser
import numpy as np

def Option_Parser(argv):

    usage='usage: %prog [options] arg\n'
    usage+='The script is to transfor scale factor information from JSON format to TH2D in ROOT file !!'
    parser = OptionParser(usage=usage)

    parser.add_option('--SF',
            type='str', dest='infile_sf', default='',
            help='Input your scale factor json file'
            )
    parser.add_option('--MC',
            type='str', dest='infile_mc', default='',
            help='Input your MC efficiency json file'
            )
    parser.add_option('--Data',
            type='str', dest='infile_data', default='',
            help='Input your data efficiency json file'
            )
    parser.add_option('-o', '--output',
            type='str', dest='output', default='output.root',
            help='Output ROOT file'
            )

    (options, args) = parser.parse_args(argv)
    return options

def abseta_pt_extractor(IDtype_content, plotname):

    for xyname in IDtype_content:
        absetaname = xyname[:xyname.find('_')]
        ptname = xyname[xyname.find('_')+1:]

        abseta_list_tmp = list()
        pt_list_tmp = list()
        value_list = list()
        error_list = list()
        for abseta_pt in IDtype_content[xyname]:
            abseta_list_tmp.append(float(abseta_pt[abseta_pt.find('[')+1:abseta_pt.find(']')].split(',')[0]))
            abseta_list_tmp.append(float(abseta_pt[abseta_pt.find('[')+1:abseta_pt.find(']')].split(',')[1]))
            for pt in IDtype_content[xyname][abseta_pt]:
                pt_list_tmp.append(float(pt[pt.find('[')+1:pt.find(']')].split(',')[0]))
                pt_list_tmp.append(float(pt[pt.find('[')+1:pt.find(']')].split(',')[1]))
        abseta_list = sorted(list(set(abseta_list_tmp)))
        pt_list = sorted(list(set(pt_list_tmp)))
        abseta_array = np.array(abseta_list)
        pt_array = np.array(pt_list)
        hist = ROOT.TH2D(plotname, plotname, len(abseta_array)-1, abseta_array, len(pt_array)-1, pt_array)

        for abseta_it in range(len(abseta_list) - 1):
            abseta_down = '%.2f' % abseta_list[abseta_it]
            abseta_up = '%.2f' % abseta_list[abseta_it+1]
            abseta_content = IDtype_content[xyname][absetaname + ':[' + abseta_down  + ',' + abseta_up + ']']

            for pt_it in range(len(pt_list) - 1):
                pt_down = '%.2f' % pt_list[pt_it]
                pt_up = '%.2f' % pt_list[pt_it+1]
                pt_content = abseta_content[ptname + ':[' + pt_down  + ',' + pt_up + ']']
                hist.SetBinContent(hist.GetBin(abseta_it+1, pt_it+1), pt_content['value'])
                hist.SetBinError(hist.GetBin(abseta_it+1, pt_it+1), pt_content['error'])
        return  hist

def hist_art_maker(raw_hist):
    raw_hist.SetOption('COLZ TEXTE')
    raw_hist.GetXaxis().SetTitle('|#eta|')
    raw_hist.GetYaxis().SetTitle('Pt (GeV)')
    raw_hist.SetTitleFont(62,"xyz");
    raw_hist.SetLabelFont(62,"xyz");
    raw_hist.SetLabelSize(0.04,"xyz");
    raw_hist.SetTitleSize(0.04,"xyz");

def json_reader(infilename):
    with open(infilename, 'r') as infile:
        text = json.load(infile)
        return text

def make_rootfile(argv):

    option = Option_Parser(argv)
    if option.infile_sf == '':
        print '[Error] : Please input your SF json file !'
        return 1
    text = json_reader(option.infile_sf)

    outfile = ROOT.TFile(option.output, 'recreate')
    ROOT.gStyle.SetOptStat(0)
    for IDtype_key in text:
        print '[INFO] : Writting %s ' % IDtype_key
        outfile.mkdir(IDtype_key)
        outfile.cd(IDtype_key)
        hist = abseta_pt_extractor(text[IDtype_key], 'abseta_pt_scalefactor')
        hist_art_maker(hist)
        hist.Write()
        if option.infile_mc != '':
            outfile.mkdir(IDtype_key + '/' + 'efficiencyMC')
            outfile.cd(IDtype_key + '/' + 'efficiencyMC')
            text_mc = json_reader(option.infile_mc)
            hist_mc = abseta_pt_extractor(text_mc[IDtype_key], 'abseta_pt_MC')
            hist_art_maker(hist_mc)
            hist_mc.Write()

        if option.infile_data != '':
            outfile.mkdir(IDtype_key + '/' + 'efficiencyData')
            outfile.cd(IDtype_key + '/' + 'efficiencyData')
            text_data = json_reader(option.infile_data)
            hist_data = abseta_pt_extractor(text_data[IDtype_key], 'abseta_pt_Data')
            hist_art_maker(hist_data)
            hist_data.Write()
    outfile.Close()
    print '[INFO] : %s file have been produced' % option.output

if __name__ == '__main__':
    sys.exit(make_rootfile(sys.argv[1:]))
