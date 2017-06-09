#!/usr/bin/python
import os, sys, math, datetime
from ROOT import gROOT, gStyle, TFile, TTree, TH1F, TH1D, TCanvas, TPad, TMath, TF1, TLegend, gPad, gDirectory, TLine, TArrow
from ROOT import kRed, kBlue, kGreen, kWhite

from datetime import datetime, date, time

from collections import OrderedDict
import numpy

sys.path.append(os.path.abspath(os.path.curdir))

from Plotter import parseInputArgs
options = parseInputArgs()

gROOT.SetBatch()
gROOT.LoadMacro("Plotter/UMDStyle.C")
from ROOT import SetUMDStyle
SetUMDStyle()

gROOT.LoadMacro("CMS_lumi.C")
#gROOT.ProcessLine(".L CMS_lumi.C+")
from ROOT import CMS_lumi, AddressOf

####################################################################################################
####################################################################################################
if __name__ == '__main__':
    myfile = {}
    myhist = {}
    myhComb = {}

    iPeriod = 0 # 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
    iPos    = 12
    writeExtraText = True

    plots_ = [
        "hmedtheta2D","hlogmedtheta2D","hmedipXYSig","hlogmedipXYSig",
        "H_T","hntrk",
        "hmass",
        "hnjet",
        "hjpta","hetaa","h_nemg"
        ]

    print "# of plots = ",len(plots_)

    labels = {}
    for iplot in plots_:
        labels[iplot] = ""

    ## specify labels
    labels["H_T"] = "H_{T} [GeV]"
    labels["hmass"] = "M [GeV]"
    for i in range(1,5):
        labels["H_T%i"%i] = "H_{T} [GeV]"
        labels["hTrig%id"%i] = "H_{T} [GeV]"
        labels["hTrig%in"%i] = "H_{T} [GeV]"

    labels["hntrk"] = "N_{trk} of jets"
    labels["h_nemg"] = "N_{tagged}"
    labels["hnjet"] = "N_{good jet}"
    labels["hetaa"] = "#eta"
    labels["hjpta"] = "Jet p_{T} [GeV]"
    labels["halpha"] = "Jet #alpha_{max}"
    labels["hacut_alphamax"] = "Jet #alpha_{max}"
    labels["hbcut_alphamax"] = "Jet #alpha_{max}"

    labels["hipXY"] = "IP^{2D} [cm]"
    labels["hmaxipXY"] = "IP^{2D}_{max} [cm]"
    labels["hipXYSig"] = "IP^{2D}_{sig}"

    labels["hmeanipXY"] = "#bar{IP}^{2D} [cm]"
    labels["hmedipXY"] = "#hat{IP}^{2D} [cm]"

    labels["hmedipXYSig"] = "#hat{IP}^{2D}_{sig}"
    labels["hmeanipXYSig"] = "#bar{IP}^{2D}_{sig}"
    labels["hmedtheta2D"] = "#hat{#Theta}_{2D}"

    labels["hlogmedipXYSig"] = "log_{10}(#hat{IP}^{2D}_{sig})"
    labels["hlogmeanipXYSig"] = "log_{10}(#bar{IP}^{2D}_{sig})"
    labels["hlogmedtheta2D"] = "log_{10}(#hat{#Theta}_{2D})"



    runtype = 1
    labeltype = 0 # 1: Val vs. Ref; 2: ntag1 vs ntag2
    norm = False

    anaTag1 = "analysis_20170523_v0_p20170608_UMD_test2p1"
    anaTag2 = "analysis_20170523_v0_p20170608_UMD_test1"
    RunLists = ["QCD80_HT1000to1500","QCD80_HT1500to2000","QCD80_HT1000to2000","QCD80_HT1000toInf","QCD80_HT1500toInf","QCD80_HT2000toInf"]
    #RunLists = ["ModelA"]

    for fileTag in RunLists:
        dirTag = fileTag.replace("_","")
        print "\n>>>>>>>> Processing: %s"%(fileTag)

        myfile["QCD"] = TFile("histos/%s/SumHists%s.root"%(anaTag1,fileTag))

        ## For runtype != 1
        myfile["QCDRef"] = TFile("histos/%s/SumHists%s.root"%(anaTag1,fileTag))
        myfile["QCDVal"] = TFile("histos/%s/SumHists%s.root"%(anaTag2,fileTag))

        colors_={}
        colors_["QCD"] = 1

        if runtype == 1:
            fTag = "%s_%s"%(anaTag1,dirTag)
        elif runtype == 2:
            fTag = "%s_vs_%s_%s_SR"%(anaTag1,anaTag2,dirTag)
        else:
            fTag = "%s_vs_%s_%s_FR"%(anaTag1,anaTag2,dirTag)

        outDir = "/home/jengbou/workspace/Data/EMJPlots/%s_%s"%(options.outtag,fTag)
        #outDir = "/Users/jengbou/hepcmsumd/workspace/CMSSW_7_6_3/src/EmergingJetAnalysis/histsQCD/EMJPlots/%s_%s"%(options.outtag,fTag)

        try:
            os.makedirs(outDir)
        except:
            pass


        for iplot in plots_:
            rebin_ = 20
            if iplot.find("H_T")!=-1:
                rebin_ = 4
                if runtype == 1:
                    myhist["QCDRef_%s"%iplot] = myfile["QCD"].Get(iplot+"3")
                    myhist["QCDVal_%s"%iplot] = myfile["QCD"].Get(iplot+"FR")
                elif runtype == 2:
                    myhist["QCDRef_%s"%iplot] = myfile["QCDRef"].Get(iplot+"3")
                    myhist["QCDVal_%s"%iplot] = myfile["QCDVal"].Get(iplot+"3")
                else:
                    myhist["QCDRef_%s"%iplot] = myfile["QCDRef"].Get(iplot+"FR")
                    myhist["QCDVal_%s"%iplot] = myfile["QCDVal"].Get(iplot+"FR")
            elif iplot.find("hmass")!=-1:
                if runtype == 1:
                    myhist["QCDRef_%s"%iplot] = myfile["QCD"].Get(iplot)
                    myhist["QCDVal_%s"%iplot] = myfile["QCD"].Get(iplot+"FR")
                elif runtype == 2:
                    myhist["QCDRef_%s"%iplot] = myfile["QCDRef"].Get(iplot)
                    myhist["QCDVal_%s"%iplot] = myfile["QCDVal"].Get(iplot)
                else:
                    myhist["QCDRef_%s"%iplot] = myfile["QCDRef"].Get(iplot+"FR")
                    myhist["QCDVal_%s"%iplot] = myfile["QCDVal"].Get(iplot+"FR")
            else:
                if iplot in ["hmedtheta2D","hjpta","hetaa","h_nemg"]: rebin_ = 1
                if iplot in ["hnjet"]: rebin_ = 1
                if iplot in ["hntrk"]: rebin_ = 5
                if iplot in ["hetaa"]: rebin_ = 4
                if iplot in ["hjpta"]: rebin_ = 2
                if iplot in ["hmedipXYSig"]: rebin_ = 5
                if runtype == 1:
                    myhist["QCDRef_%s"%iplot] = myfile["QCD"].Get(iplot+"SR")
                    myhist["QCDVal_%s"%iplot] = myfile["QCD"].Get(iplot+"FR")
                elif runtype == 2:
                    myhist["QCDRef_%s"%iplot] = myfile["QCDRef"].Get(iplot+"SR")
                    myhist["QCDVal_%s"%iplot] = myfile["QCDVal"].Get(iplot+"SR")
                else:
                    myhist["QCDRef_%s"%iplot] = myfile["QCDRef"].Get(iplot+"FR")
                    myhist["QCDVal_%s"%iplot] = myfile["QCDVal"].Get(iplot+"FR")

            cvsName = "Fig_%s"%(iplot)
            fnameTag = "%s/%s"%(outDir,cvsName)
            c1 = TCanvas(fnameTag,fnameTag,800,800)
            pad1 = TPad("pad1", "", 0, 0.3, 1, 1)
            pad2 = TPad("pad2", "", 0, 0, 1, 0.30)
            pad1.SetBottomMargin(0.02)
            #pad2.SetTopMargin(0)
            pad2.SetBottomMargin(0.4)
            pad2.SetGridy(1)
            pad1.Draw()
            pad2.Draw()

            pad1.cd()

            if myhist["QCDRef_%s"%iplot] == None: continue

            print "# events ref = ", myhist["QCDRef_%s"%iplot].Integral(0,myhist["QCDRef_%s"%iplot].GetNbinsX()+2)
            print "# events val = ", myhist["QCDVal_%s"%iplot].Integral(0,myhist["QCDVal_%s"%iplot].GetNbinsX()+2)
            #myhist["QCD_%s"%iplot].Rebin(rebin_)
            myhist["QCDRef_%s"%iplot].Rebin(rebin_)
            myhist["QCDVal_%s"%iplot].Rebin(rebin_)

    ##        myhist["QCDRef_%s"%iplot].SetLineColor(2)
    ##        myhist["QCDRef_%s"%iplot].Draw("hist")

    ##        myhist["QCDRef_%s"%iplot].GetXaxis().SetRangeUser(0.,1.)
    ##        myhist["QCDRef_%s"%iplot].GetYaxis().SetRangeUser(1.e-3,5.)
    ##        gPad.Update()

            myhist["QCDRef_%s"%iplot].GetXaxis().SetLabelOffset(0.02)
            myhist["QCDRef_%s"%iplot].SetLineColor(2)
            myhist["QCDRef_%s"%iplot].SetMarkerColor(2)
            myhist["QCDRef_%s"%iplot].SetMarkerSize(0)

            if (norm and iplot in ["hnjet","hetaa","hjpta","hntrk","hmass"]) or labeltype == 2:
                myhist["QCDRef_%s"%iplot].Scale(1./myhist["QCDRef_%s"%iplot].Integral())
                myhist["QCDRef_%s"%iplot].GetYaxis().SetTitle("A.U.")
            else:
                myhist["QCDRef_%s"%iplot].GetYaxis().SetTitle("Events")

            myhist["QCDRef_%s"%iplot].Draw("hist e")

            myhist["QCDVal_%s"%iplot].SetLineColor(4)
            myhist["QCDVal_%s"%iplot].SetMarkerColor(4)
            myhist["QCDVal_%s"%iplot].SetMarkerSize(0)
            if (norm and iplot in ["hnjet","hetaa","hjpta","hntrk","hmass"]) or labeltype == 2:
                myhist["QCDVal_%s"%iplot].Scale(1./myhist["QCDVal_%s"%iplot].Integral())

            myhist["QCDVal_%s"%iplot].Draw("hist e sames")


            if iplot in ["H_T"]:
                xMin_ = 1000.0
                xMax_ = 3000.0
                myhist["QCDRef_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                myhist["QCDVal_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                gPad.Update()
                c1.Update()

            if iplot in ["hmedtheta2D"]:
                xMin_ = 0.0
                xMax_ = 0.05
                myhist["QCDRef_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                myhist["QCDVal_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                gPad.Update()
                c1.Update()

            if iplot in ["hmedipXYSig"]:
                xMin_ = 0.0
                xMax_ = 10.0
                myhist["QCDRef_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                myhist["QCDVal_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                gPad.Update()
                c1.Update()
    ##        myhist["QCDVal_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 1.1*valYmax)

            if iplot in ["hlogmedipXYSig"]:
                xMin_ = -1.0
                xMax_ = 2.0
                myhist["QCDRef_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                myhist["QCDVal_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                gPad.Update()
                c1.Update()

            if iplot in ["hlogmedtheta2D"]:
                xMin_ = -3.0
                xMax_ = -1.0
                myhist["QCDRef_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                myhist["QCDVal_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                gPad.Update()
                c1.Update()


            if iplot in ["h_nemg"]:
                xMin_ = 0.0
                xMax_ = 4.0
                myhist["QCDRef_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                myhist["QCDVal_%s"%iplot].GetXaxis().SetRangeUser(xMin_, xMax_)
                gPad.Update()
                c1.Update()

            ###############################################################
            valYmax =  -99999.
            valYmin = 99999.
            logY_ = False

            for nh in [myhist["QCDRef_%s"%iplot],myhist["QCDVal_%s"%iplot]]:
                if iplot in ["h_nemg","hnjet"]: logY_ = True
                for i in range(1,nh.GetNbinsX()+1):
                    if  nh.GetBinContent(i) < valYmin: valYmin = nh.GetBinContent(i)
                    if valYmin <= 0. or logY_:
                        if iplot in ["h_nemg"] and labeltype != 2:
                            valYmin = 1.
                        else:
                            valYmin = 1.e-4
                    if nh.GetBinContent(i) > valYmax:
                        valYmax = nh.GetBinContent(i)
            #print valYmax
            if not logY_:
                myhist["QCDRef_%s"%iplot].GetYaxis().SetRangeUser(0., 1.2*valYmax)
                myhist["QCDVal_%s"%iplot].GetYaxis().SetRangeUser(0., 1.2*valYmax)
            else:
                myhist["QCDRef_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 3.*valYmax)
                myhist["QCDVal_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 3.*valYmax)

            gPad.Update()
            c1.Update()


            ###############################################################
            leg = TLegend(0.72,0.7,0.94,0.91)
            leg.SetFillColor(kWhite)
            leg.SetLineColor(kWhite)

            if runtype == 1:
                leg.AddEntry(myhist["QCDRef_%s"%iplot],"QCD MC","l")
                leg.AddEntry(myhist["QCDVal_%s"%iplot],"Fake bkg","l")
            elif runtype == 2 and labeltype == 2:
                leg.AddEntry(myhist["QCDRef_%s"%iplot],"N_{tag}=1 (SR)","l")
                leg.AddEntry(myhist["QCDVal_%s"%iplot],"N_{tag}=2 (SR)","l")
            elif runtype == 3 and labeltype == 2:
                leg.AddEntry(myhist["QCDRef_%s"%iplot],"N_{tag}=1 (FR)","l")
                leg.AddEntry(myhist["QCDVal_%s"%iplot],"N_{tag}=2 (FR)","l")
            else:
                leg.AddEntry(myhist["QCDRef_%s"%iplot],"Ref","l")
                leg.AddEntry(myhist["QCDVal_%s"%iplot],"Val","l")

            leg.Draw()

            if logY_:
                gPad.SetLogy()
                gPad.Update()
                c1.Update()


            ## ratio plot:
            pad2.cd()
            ## FIXME
            htmp = myhist["QCDRef_%s"%iplot].Clone("ratio_den_%s"%iplot)
            rtemp = myhist["QCDVal_%s"%(iplot)].Clone("ratio_num_%s"%(iplot))
            #if iplot.find("alpha")!=-1: rtemp.GetXaxis().SetRangeUser(0.,1.)
            rtemp.GetXaxis().SetTitle(labels[iplot])
            #rtemp = myhist["QCD_%s"%iplot].Clone("ratio_%s"%(iplot))
            rtemp.Divide(htmp)
            rtemp.SetMarkerStyle(1)
            rtemp.SetMarkerColor(4)
            rtemp.GetXaxis().SetTitleSize(0.12)
            rtemp.GetYaxis().SetTitleSize(0.12)
            rtemp.GetXaxis().SetTitleOffset(1.20)
            rtemp.GetYaxis().SetTitleOffset(0.60)
            rtemp.GetXaxis().SetLabelSize(0.1)
            rtemp.GetYaxis().SetLabelSize(0.1)
            rtemp.GetYaxis().SetNdivisions(5)
            if labeltype == 1:
                rtemp.GetYaxis().SetTitle("#frac{Val}{Ref}")
            elif labeltype == 2:
                rtemp.GetYaxis().SetTitle("#frac{N_{tag}=2}{N_{tag}=1}")
            else:
                rtemp.GetYaxis().SetTitle("#frac{FakeBkg}{MC}")
            rtemp.GetYaxis().SetRangeUser(0.,2.)
            rtemp.DrawCopy("c")
            del htmp
            if runtype == 1 and fileTag == "QCD80_HT1000toInf" and iplot == "hnjet":
                for ibin in range(rtemp.GetNbinsX()+1):
                    if rtemp.GetBinContent(ibin+1)>0:
                        #print "Ratio bin [%i] = %7.5f"%(ibin+1,1./rtemp.GetBinContent(ibin+1))
                        print "Ratio [%5.2f] = %7.5f"%(rtemp.GetBinLowEdge(ibin+1),1./rtemp.GetBinContent(ibin+1))
                    else:
                        #print "Ratio bin [%i] = %7.5f"%(ibin+1,1.)
                        print "Ratio [%5.2f] = %7.5f"%(rtemp.GetBinLowEdge(ibin+1),1.)

            CMS_lumi(c1,iPeriod,iPos);
            gPad.RedrawAxis()
            c1.SaveAs("%s.pdf"%fnameTag)
            c1.SaveAs("%s.png"%fnameTag)


