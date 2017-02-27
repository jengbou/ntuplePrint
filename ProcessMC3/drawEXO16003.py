#!/usr/bin/python
import os, sys, math, datetime
from ROOT import gROOT, gStyle, TFile, TTree, TH1F, TH1D, TCanvas, TPad, TMath, TF1, TLegend, gPad, gDirectory, TLine, TArrow
from ROOT import kRed, kBlue, kGreen, kWhite

from datetime import datetime, date, time

from collections import OrderedDict
import numpy

sys.path.append(os.path.abspath(os.path.curdir))

from Plotter import parseLYAnaInputArgs
options = parseLYAnaInputArgs()

from Plotter.CommonTools import CalcD, AlphaSourceFitter, CalcNPE, HistToTGRaphErrorsUniform

gROOT.SetBatch()
gROOT.LoadMacro("Plotter/AtlasStyle.C")
from ROOT import SetAtlasStyle
SetAtlasStyle()

gROOT.LoadMacro("CMS_lumi.C")
#gROOT.ProcessLine(".L CMS_lumi.C+")
from ROOT import CMS_lumi, AddressOf

####################################################################################################
####################################################################################################
if __name__ == '__main__':
    myfile  = {}
    myhist  = {}
    myhComb = {}

    iPeriod = 0 # 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
    iPos    = 12
    writeExtraText = True

    plots_ = [
        "count","acount","h_nemg",
##        "hjetcut","hjetchf",
##        "hnjet","hpt","heta","heta2",
        "H_T","H_T1","H_T2","H_T3","H_T4",
##        "hpt1","hpt2","hpt3",
##        "hpt4","hbcut_ntrkpt1","hacut_ntrkpt1","hbcut_nef","hacut_nef",
##        "hbcut_cef","hacut_cef","hbcut_alphamax","hacut_alphamax","hHTnm1",
##        "hpt1nm1","hpt2nm1","hpt3nm1","hpt4nm1","halphanm1",
##        "hmaxipnm1","hnHitsnm1","hntrk1nm1","hnemnm1",
##        "hnmaxipnm1","hn2maxipnm1","hjptfrb","hjptfra1",
##        "hjptfra2","hjptfrbc","hjptfra1c","hjptfra2c","hjptb",
##        "hjpta","haMgj","hHTko","hpt1ko","hpt2ko",
##        "hpt3ko","hpt4ko",
        "hmass","halpha","halphaPS",
        "hTrig1d","hTrig1n","hTrig2d","hTrig2n","hTrig3d","hTrig3n","h_ntag",
        "hmedtheta2DPS","hlogmedtheta2DPS","hmedipXYSigPS","hlogmedipXYSigPS",
        "hmedtheta2DSR","hlogmedtheta2DSR","hmedipXYSigSR","hlogmedipXYSigSR",
        ]

    print "# of plots = ",len(plots_)

    plotsComb_ = [
        ["hipXYEJ","hipXYnEJ"],
        ["hipXYSigEJ","hipXYSignEJ"],
        ["hmaxipXYEJ","hmaxipXYnEJ"],
        ["hmeanipXYEJ","hmeanipXYnEJ"],
        ["hlogmedipXYSigEJ","hlogmedipXYSignEJ"],
        ["hlogmeanipXYSigEJ","hlogmeanipXYSignEJ"],
        ["hmedipXYSigEJ","hmedipXYSignEJ"],
        ["hmeanipXYSigEJ","hmeanipXYSignEJ"],
        ["hmedipXYEJ","hmedipXYnEJ"],
        ["hmedtheta2DEJ","hmedtheta2DnEJ"],
        ["hlogmedtheta2DEJ","hlogmedtheta2DnEJ"],
        ]

    print "# of plots (comb.x2) = ",len(plotsComb_*2)

    ## Plots to be shown in log scale
    ## match plot with pattern
    logYplots = ["alpha","hipXY","hmeanipXY","hmedipXY","hmeanipXYSig","hmedipXYSig","hmaxipXY","hipXYSig","hmedtheta2D"]
    ## match plot with exact name:
    logYplotsE = ["H_T1","H_T2","H_T3","H_T4"]

    labels = {}
    for iplot in plots_:
        labels[iplot] = ""
    for iplot in plotsComb_:
        labels[iplot[0]] = ""
        labels[iplot[1]] = ""

    ## specify labels
    labels["H_T"] = "H_{T} [GeV]"
    for i in range(1,5):
        labels["H_T%i"%i] = "H_{T} [GeV]"
        labels["hTrig%id"%i] = "H_{T} [GeV]"
        labels["hTrig%in"%i] = "H_{T} [GeV]"

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

    myfile["ModelA"] = TFile("histos/analysis_20170223_v0_p20170226/SumHistsModelA.root")
    myfile["ModelB"] = TFile("histos/analysis_20170223_v0_p20170226/SumHistsModelB.root")
    myfile["QCD"] = TFile("histos/analysis_20170223_v0_p20170226/SumHistsQCD74.root")

    colors_={}
    colors_["QCD"] = 1
    colors_["ModelA"] = 2
    colors_["ModelB"] = 4

    #fTag = datetime.now().strftime("%Y%m%d_%H%M%S")
    fTag = "analysis20170223_v0" #"analysis20170215_v0_p20170222"
    #fTag = "20170216_203317_v0p1"
    outDir = "/data/users/jengbou/workspace/Data/EMJPlots/%s_%s"%(options.outtag,fTag)

    try:
        os.makedirs(outDir)
    except:
        pass

    ####################################################################
    # Normal plots
    ####################################################################
    for iplot in plots_:
        hNameTmp = iplot.rstrip("PS").rstrip("SR")
        print "hNameTmp = ",hNameTmp
        rebin_ = 1
        updateXrange = False
        updateYrange = False
        logY_ = False
        if hNameTmp in [
            "H_T","H_T1","H_T2","H_T3","H_T4",
            "hmeanipXYSig","hmedipXYSig","halpha","hmedtheta2D",
            "hTrig1d","hTrig1n","hTrig2d","hTrig2n","hTrig3d","hTrig3n"
            ]:
            rebin_ = 1
        elif hNameTmp in ["hipXYSig","hipXY"]:
            rebin_ = 2
        else:
            rebin_ = 10

        ## change only SR
        if iplot in ["hmedipXYSigSR"]:
            rebin_ = 5

        if hNameTmp in logYplots+logYplotsE: logY_ = True

        myhist["ModelA_%s"%iplot] = myfile["ModelA"].Get(iplot)
        myhist["ModelB_%s"%iplot] = myfile["ModelB"].Get(iplot)
        myhist["QCD_%s"%iplot] = myfile["QCD"].Get(iplot)

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

        if iplot.find("count")==-1:
            myhist["QCD_%s"%iplot].Rebin(rebin_)
            myhist["ModelA_%s"%iplot].Rebin(rebin_)
            myhist["ModelB_%s"%iplot].Rebin(rebin_)

        ## Model A
        if myhist["ModelA_%s"%iplot].Integral()!=0 and iplot.find("count")==-1:
            myhist["ModelA_%s"%iplot].Scale(1./myhist["ModelA_%s"%iplot].Integral())
        myhist["ModelA_%s"%iplot].SetLineColor(2)
        myhist["ModelA_%s"%iplot].GetXaxis().SetTitle(labels[hNameTmp])
        myhist["ModelA_%s"%iplot].GetYaxis().SetTitle("A.U.")
        myhist["ModelA_%s"%iplot].GetXaxis().SetLabelOffset(0.02)
        myhist["ModelA_%s"%iplot].Draw()

        ## Model B
        if myhist["ModelB_%s"%iplot].Integral()!=0 and iplot.find("count")==-1:
            myhist["ModelB_%s"%iplot].Scale(1./myhist["ModelB_%s"%iplot].Integral())
        myhist["ModelB_%s"%iplot].SetLineColor(4)
        myhist["ModelB_%s"%iplot].Draw("sames")

        ## QCD
        if myhist["QCD_%s"%iplot].Integral()!=0 and iplot.find("count")==-1:
            myhist["QCD_%s"%iplot].Scale(1./myhist["QCD_%s"%iplot].Integral())
        myhist["QCD_%s"%iplot].GetXaxis().SetLabelOffset(0.02)
        myhist["QCD_%s"%iplot].Draw("sames")


        leg = TLegend(0.72,0.7,0.94,0.91)
        leg.SetFillColor(kWhite)
        leg.SetLineColor(kWhite)

        leg.AddEntry(myhist["QCD_%s"%iplot],"QCD","l")
        leg.AddEntry(myhist["ModelA_%s"%iplot],"Model A","l")
        leg.AddEntry(myhist["ModelB_%s"%iplot],"Model B","l")
        leg.Draw()


        ###############################################################
        # Calculate better y range to draw plots
        ###############################################################
        valYmax =  -9999.
        valYmin = 9999.

        for nh in ["ModelA_%s"%(iplot),"ModelB_%s"%(iplot),"QCD_%s"%(iplot)]:
            if (nh not in myhist.keys()):
                print "%s not found"%nh
                continue
            for i in range(1,myhist[nh].GetNbinsX()+1):
                if  myhist[nh].GetBinContent(i) < valYmin: valYmin = 0.9*myhist[nh].GetBinContent(i)
                if valYmin <= 0. and logY_:
                    valYmin = 1.e-4
                if valYmin < 0: valYmin = 1.e-6
                if myhist[nh].GetBinContent(i) > valYmax:
                    valYmax = myhist[nh].GetBinContent(i)


        ###############################################################
        if iplot in ["hTrig1d","hTrig1n","hTrig2d","hTrig2n","hTrig3d","hTrig3n"]: valYmin = 0.
        myhist["QCD_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 1.2*valYmax)
        myhist["ModelA_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 1.2*valYmax)
        myhist["ModelB_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 1.2*valYmax)
        #print "Ymin,Ymax = %f,%f"%(valYmin,1.2*valYmax)
        gPad.Update()

        ###############################################################
        ## Change plot range for specified plots
        ## Make sure to cross check overlap of hNameTmp and iplot items
        ###############################################################
        ## find plots with exact name
        if iplot in ["H_T","H_T1","H_T2","H_T3","H_T4"]:
            updateXrange = True
            xMin_ = 0.0
            xMax_ = 3500.0
            if iplot in ["H_T"]:
                valYmin = 1.e-6
                valYmax *= 10
            else:
                valYmin = 1.e-4
                valYmax *= 1.5

        if iplot in ["hmedipXYSigPS"]:
            updateXrange = True
            xMin_ = 0.0
            xMax_ = 10.0

        if iplot in ["hmedipXYSigSR"]:
            updateXrange = True
            xMin_ = 0.0
            xMax_ = 100.0

        ## find plots with pattern ==> make sure to cross check any overlap w.r.t. items specified by exact name
        if hNameTmp in ["hipXY","hmeanipXY","hmedipXY","hmeanipXYSig","hmedipXYSig","hmaxipXY","hipXYSig","hmedtheta2D"]:
            updateYrange = True
            if valYmin < 1.e-4: valYmin = 1.e-4
            if hNameTmp in ["hmeanipXY"]: valYmin = 1.e-4
            if hNameTmp in ["hmedipXY"]: valYmin = 1.e-5
            valYmax *= 1.5;
            if hNameTmp in ["hmedtheta2D"]:
                updateXrange = True
                xMin_ = 0.0
                xMax_ = 0.2

            elif hNameTmp in ["hipXY"]:
                xMin_ = 0.0
                xMax_ = 1.0
                updateXrange = True

            elif hNameTmp in ["hipXYSig","hmeanipXYSig"]:
                xMin_ = 0.0
                xMax_ = 10.0
                updateXrange = True

        ###############################################################
        # Add cut line for jet tagger plots
        ###############################################################
        if iplot.find("alpha")!=-1:
            myhist["QCD_%s"%iplot].GetXaxis().SetRangeUser(0.,1.)
            myhist["ModelA_%s"%iplot].GetXaxis().SetRangeUser(0.,1.)
            myhist["ModelB_%s"%iplot].GetXaxis().SetRangeUser(0.,1.)

            myhist["QCD_%s"%iplot].GetYaxis().SetRangeUser(1.e-3,5.)
            myhist["ModelA_%s"%iplot].GetYaxis().SetRangeUser(1.e-3,5.)
            myhist["ModelB_%s"%iplot].GetYaxis().SetRangeUser(1.e-3,5.)

            gPad.Update()

            reflcut = TLine()
            reflcut.SetLineStyle(2)
            reflcut.SetLineWidth(2)
            reflcut.SetLineColor(9)
            reflcut.DrawLine(0.05,valYmin,0.05,1.5*valYmax)

            refacut = TArrow(0.05,1.15*valYmax,0.015,1.15*valYmax,0.005)
            refacut.SetLineWidth(2)
            refacut.SetLineColor(9)
            refacut.Draw()

            if valYmin < 1.e-3: valYmin = 1.e-3
            if 1.5*valYmax < 5.0: valYmax = 5.0
            updateYrange = True
            logY_ = True


        ###############################################################
        # Apply ranges update
        ###############################################################

        if updateXrange:
            myhist["QCD_%s"%(iplot)].GetXaxis().SetRangeUser(xMin_,xMax_)
            myhist["ModelA_%s"%(iplot)].GetXaxis().SetRangeUser(xMin_,xMax_)
            myhist["ModelB_%s"%(iplot)].GetXaxis().SetRangeUser(xMin_,xMax_)

        if updateYrange:
            myhist["QCD_%s"%(iplot)].GetYaxis().SetRangeUser(valYmin,valYmax)
            myhist["ModelA_%s"%(iplot)].GetYaxis().SetRangeUser(valYmin,valYmax)
            myhist["ModelB_%s"%(iplot)].GetYaxis().SetRangeUser(valYmin,valYmax)

        if (updateXrange or updateYrange):
            gPad.Update()
            c1.Update()

        if logY_:
            gPad.SetLogy()
            gPad.Update()
            c1.Update()


        ############################
        # This is special case
        ############################
        if iplot.find("count")!=-1:
            myhist["QCD_%s"%iplot].Draw("hist text0")
            myhist["ModelA_%s"%iplot].Draw("sames hist text0")
            myhist["ModelB_%s"%iplot].Draw("sames hist text0")
            leg.Draw()
            myhist["QCD_%s"%iplot].GetXaxis().SetRangeUser(0.,10.)
            myhist["ModelA_%s"%iplot].GetXaxis().SetRangeUser(0.,10.)
            myhist["ModelB_%s"%iplot].GetXaxis().SetRangeUser(0.,10.)

            gPad.Update()

            if valYmin < 1.e-4 :valYmin = 1.e-4
##            if iplot=="acount":
            myhist["QCD_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 300.*valYmax)
            myhist["ModelA_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 300.*valYmax)
            myhist["ModelB_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 300.*valYmax)
            #print "Ymin,Ymax = %f,%f"%(valYmin,1.5*valYmax)
            gPad.SetLogy()
##            else:
##                myhist["QCD_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 1.6*valYmax)
##                myhist["ModelA_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 1.6*valYmax)
##                myhist["ModelB_%s"%iplot].GetYaxis().SetRangeUser(valYmin, 1.6*valYmax)

            gPad.Update()
            c1.Update()


        ############################
        # Ratio plots
        ############################
        pad2.cd()
        ## FIXME
        htmp = myhist["ModelA_%s"%iplot].Clone("ratio_den_%s"%iplot)
        rtemp = myhist["ModelA_%s"%(iplot)].Clone("ratio_num_%s"%(iplot))
        for i in xrange(rtemp.GetNbinsX()):
            rtemp.SetBinContent(i+1,1)
        #rtemp.Divide(htmp)
        if iplot.find("alpha")!=-1: rtemp.GetXaxis().SetRangeUser(0.,1.)
        rtemp.GetXaxis().SetTitle(labels[hNameTmp])
        rtemp.GetXaxis().SetTitleSize(0.12)
        rtemp.GetYaxis().SetTitleSize(0.12)
        rtemp.GetXaxis().SetTitleOffset(1.20)
        rtemp.GetYaxis().SetTitleOffset(0.60)
        rtemp.GetXaxis().SetLabelSize(0.1)
        rtemp.GetYaxis().SetLabelSize(0.1)
        rtemp.GetYaxis().SetNdivisions(5)
        rtemp.GetXaxis().SetTitle(labels[hNameTmp])
        rtemp.GetYaxis().SetTitle("#frac{Data-MC}{MC}")
        rtemp.GetYaxis().SetRangeUser(0.,2.)
        rtemp.DrawCopy("hist")

        CMS_lumi(c1,iPeriod,iPos);
        gPad.RedrawAxis()
        c1.SaveAs("%s.pdf"%fnameTag)
        c1.SaveAs("%s.png"%fnameTag)


    ####################################################################
    # Combined plots
    ####################################################################
    for iplot in plotsComb_:
        updateXrange = False
        updateYrange = False
        logY_ = False

        if iplot[0] in ["hmeanipXYSigEJ","hmedipXYSigEJ","hmedtheta2DEJ"]:
            rebin_ = 1
        elif iplot[0] in ["hipXYSigEJ","hipXYEJ"]:
            rebin_ = 2
        else:
            rebin_ = 10

        #print iplot[0]
        hNameTmp = iplot[0][0:iplot[0].find("EJ")]
        print "Histo Name: %s ; rebin = %i"%(hNameTmp,rebin_)
        myhComb["ModelA_%s_SR"%(hNameTmp)] = myfile["ModelA"].Get(iplot[0])
        myhComb["ModelB_%s_SR"%(hNameTmp)] = myfile["ModelB"].Get(iplot[0])
        myhComb["QCD_%s_SR"%(hNameTmp)] = myfile["QCD"].Get(iplot[0])

        myhComb["ModelA_%s_SB"%(hNameTmp)] = myfile["ModelA"].Get(iplot[1])
        myhComb["ModelB_%s_SB"%(hNameTmp)] = myfile["ModelB"].Get(iplot[1])
        myhComb["QCD_%s_SB"%(hNameTmp)] = myfile["QCD"].Get(iplot[1])

        cvsName = "Fig_Comb_%s"%(hNameTmp)
        fnameTag = "%s/%s"%(outDir,cvsName)
        c2 = TCanvas(fnameTag,fnameTag,800,800)
        pad1 = TPad("pad1", "", 0, 0.3, 1, 1)
        pad2 = TPad("pad2", "", 0, 0, 1, 0.30)
        pad1.SetBottomMargin(0.02)
        #pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.4)
        pad2.SetGridy(1)
        pad1.Draw()
        pad2.Draw()

        pad1.cd()

        myhComb["QCD_%s"%(hNameTmp)] = myhComb["QCD_%s_SR"%(hNameTmp)].Clone("QCD_%s"%(hNameTmp))
        myhComb["ModelA_%s"%(hNameTmp)] = myhComb["ModelA_%s_SR"%(hNameTmp)].Clone("ModelA_%s"%(hNameTmp))
        myhComb["ModelB_%s"%(hNameTmp)] = myhComb["ModelB_%s_SR"%(hNameTmp)].Clone("ModelB_%s"%(hNameTmp))

        myhComb["QCD_%s"%(hNameTmp)].Add(myhComb["QCD_%s_SB"%(hNameTmp)])
        myhComb["ModelA_%s"%(hNameTmp)].Add(myhComb["ModelA_%s_SB"%(hNameTmp)])
        myhComb["ModelB_%s"%(hNameTmp)].Add(myhComb["ModelB_%s_SB"%(hNameTmp)])

        myhComb["QCD_%s"%(hNameTmp)].Rebin(rebin_)
        myhComb["ModelA_%s"%(hNameTmp)].Rebin(rebin_)
        myhComb["ModelB_%s"%(hNameTmp)].Rebin(rebin_)

        if myhComb["QCD_%s"%(hNameTmp)].Integral()!=0 and hNameTmp.find("count")==-1:
            myhComb["QCD_%s"%(hNameTmp)].Scale(1./myhComb["QCD_%s"%(hNameTmp)].Integral())
            myhComb["QCD_%s"%(hNameTmp)].GetXaxis().SetTitle(labels[hNameTmp])
            myhComb["QCD_%s"%(hNameTmp)].GetYaxis().SetTitle("A.U.")
            myhComb["QCD_%s"%(hNameTmp)].GetXaxis().SetLabelOffset(0.02)
            myhComb["QCD_%s"%(hNameTmp)].Draw()

        if myhComb["ModelA_%s"%(hNameTmp)].Integral()!=0 and hNameTmp.find("count")==-1:
            myhComb["ModelA_%s"%(hNameTmp)].Scale(1./myhComb["ModelA_%s"%(hNameTmp)].Integral())
        myhComb["ModelA_%s"%(hNameTmp)].SetLineColor(2)
##        myhComb["ModelA_%s"%(hNameTmp)].GetYaxis().SetTitle("A.U.")
##        myhComb["ModelA_%s"%(hNameTmp)].GetXaxis().SetLabelOffset(0.02)
        myhComb["ModelA_%s"%(hNameTmp)].Draw("sames")

        if myhComb["ModelB_%s"%(hNameTmp)].Integral()!=0 and hNameTmp.find("count")==-1:
            myhComb["ModelB_%s"%(hNameTmp)].Scale(1./myhComb["ModelB_%s"%(hNameTmp)].Integral())
        myhComb["ModelB_%s"%(hNameTmp)].SetLineColor(4)
        myhComb["ModelB_%s"%(hNameTmp)].Draw("sames")

        ## Change plot range for specified plots
        xMin_,xMax_ = myhComb["ModelB_%s"%(hNameTmp)].GetXaxis().GetXmin(),myhComb["ModelB_%s"%(hNameTmp)].GetXaxis().GetXmax()
        updateRange = True
        if hNameTmp in ["hipXY"]:
            xMin_ = 0.0
            xMax_ = 1.0
        elif hNameTmp in ["hipXYSig","hmeanipXYSig"]:
            xMin_ = 0.0
            xMax_ = 10.0
        elif hNameTmp in ["hmedtheta2D"]:
            xMin_ = 0.0
            xMax_ = 0.1
        else: updateRange = False

        if updateRange:
            myhComb["QCD_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)
            myhComb["ModelA_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)
            myhComb["ModelB_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)
            gPad.Update()
            c2.Update()


        valYmax =  -99999.
        valYmin = 99999.

        for nh in ["ModelA_%s"%(hNameTmp),"ModelB_%s"%(hNameTmp),"QCD_%s"%(hNameTmp)]:
            if (nh not in myhComb.keys()): continue
            for i in range(1,myhComb[nh].GetNbinsX()+1):
                #print "[%s bin %i = %f]"%(nh,i,myhComb[nh].GetBinContent(i))
                if  myhComb[nh].GetBinContent(i) < valYmin: valYmin = 0.9*myhComb[nh].GetBinContent(i)
                if valYmin <= 0.:
                    valYmin = 1.e-6
                if myhComb[nh].GetBinContent(i) > valYmax:
                    valYmax = myhComb[nh].GetBinContent(i)

        #print "Ymin,Ymax = %10f,%10f"%(valYmin,1.2*valYmax)
        myhComb["QCD_%s"%(hNameTmp)].GetYaxis().SetRangeUser(valYmin, 1.2*valYmax)
        myhComb["ModelA_%s"%(hNameTmp)].GetYaxis().SetRangeUser(valYmin, 1.2*valYmax)
        myhComb["ModelB_%s"%(hNameTmp)].GetYaxis().SetRangeUser(valYmin, 1.2*valYmax)
        gPad.Update()
        c2.Update()

        leg = TLegend(0.72,0.7,0.94,0.91)
        leg.SetFillColor(kWhite)
        leg.SetLineColor(kWhite)

        leg.AddEntry(myhComb["QCD_%s"%(hNameTmp)],"QCD","l")
        leg.AddEntry(myhComb["ModelA_%s"%(hNameTmp)],"Model A","l")
        leg.AddEntry(myhComb["ModelB_%s"%(hNameTmp)],"Model B","l")
        leg.Draw()

        if hNameTmp.find("logmedipXYSig")!=-1:
            reflcut = TLine()
            reflcut.SetLineStyle(2)
            reflcut.SetLineWidth(3)
            reflcut.SetLineColor(9)
            reflcut.DrawLine(1.5,valYmin,1.5,0.8*valYmax)

            refacut = TArrow(1.51,0.5*valYmax,2.5,0.5*valYmax,0.01)
            refacut.SetLineWidth(2)
            refacut.SetLineColor(9)
            refacut.Draw()

        elif hNameTmp.find("logmedtheta2D")!=-1:
            reflcut = TLine()
            reflcut.SetLineStyle(2)
            reflcut.SetLineWidth(3)
            reflcut.SetLineColor(9)
            reflcut.DrawLine(-1.6,valYmin,-1.6,0.8*valYmax)

            refacut = TArrow(-1.49,0.65*valYmax,-0.6,0.65*valYmax,0.01)
            refacut.SetLineWidth(2)
            refacut.SetLineColor(9)
            refacut.Draw()

        gPad.RedrawAxis()

        if hNameTmp in ["hipXY","hmeanipXY","hmedipXY","hmeanipXYSig","hmedipXYSig","hmaxipXY","hipXYSig","hmedtheta2D"]:
            logY_ = True
            updateYrange = True
            valYmax *= 1.5

            if valYmin < 1.e-4: valYmin = 1.e-4
            if hNameTmp in ["hmeanipXY"]: valYmin = 1.e-4
            if hNameTmp in ["hmedipXY","hmedtheta2D"]: valYmin = 1.e-5

        if hNameTmp in ["hipXY"]:
            updateXrange = True
            xMin_ = 0.0
            xMax_ = 1.0

        elif hNameTmp in ["hipXYSig","hmeanipXYSig","hmedipXYSig"]:
            updateXrange = True
            xMin_ = 0.0
            xMax_ = 10.0

##        updateRange = True
##        if hNameTmp in ["hipXY"]:
##            xMin_ = 0.0
##            xMax_ = 1.0

##        elif hNameTmp in ["hipXYSig","hmeanipXYSig","hmedipXYSig"]:
##            xMin_ = 0.0
##            xMax_ = 10.0

##        else: UpdateRange = False

##        if updateRange:
##            myhComb["QCD_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)
##            myhComb["ModelA_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)
##            myhComb["ModelB_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)
##            gPad.Update()
##            c2.Update()

##        if logY_:
##            gPad.SetLogy()
##            gPad.Update()
##            c2.Update()

        ###############################################################
        # Apply ranges update
        ###############################################################

        if updateXrange:
            myhComb["QCD_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)
            myhComb["ModelA_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)
            myhComb["ModelB_%s"%(hNameTmp)].GetXaxis().SetRangeUser(xMin_,xMax_)

        if updateYrange:
            myhComb["QCD_%s"%hNameTmp].GetYaxis().SetRangeUser(valYmin, valYmax)
            myhComb["ModelA_%s"%hNameTmp].GetYaxis().SetRangeUser(valYmin, valYmax)
            myhComb["ModelB_%s"%hNameTmp].GetYaxis().SetRangeUser(valYmin, valYmax)

        if (updateXrange or updateYrange):
            gPad.Update()
            c2.Update()

        if logY_:
            gPad.SetLogy()
            gPad.Update()
            c2.Update()


        ############################
        # Ratio plots
        ############################
        pad2.cd()

        htmp = myhComb["ModelA_%s"%(hNameTmp)].Clone("temp_%s"%(hNameTmp))
        rtemp = myhComb["ModelA_%s"%(hNameTmp)].Clone("ratio_%s"%(hNameTmp))
        ## FIXME
        for i in xrange(rtemp.GetNbinsX()):
            rtemp.SetBinContent(i+1,1)
        #rtemp.Divide(htmp)
        rtemp.GetXaxis().SetTitleSize(0.12)
        rtemp.GetYaxis().SetTitleSize(0.12)
        rtemp.GetXaxis().SetTitleOffset(1.20)
        rtemp.GetYaxis().SetTitleOffset(0.60)
        rtemp.GetXaxis().SetLabelSize(0.1)
        rtemp.GetYaxis().SetLabelSize(0.1)
        rtemp.GetYaxis().SetNdivisions(5)
        rtemp.GetXaxis().SetTitle(labels[hNameTmp])
        rtemp.GetYaxis().SetTitle("#frac{Data-MC}{MC}")
        rtemp.GetYaxis().SetRangeUser(0.,2.)

        if updateRange: rtemp.GetXaxis().SetRangeUser(xMin_,xMax_)

        rtemp.DrawCopy("c")


        CMS_lumi(c2,iPeriod,iPos);
        gPad.RedrawAxis()
        c2.SaveAs("%s.pdf"%fnameTag)
        c2.SaveAs("%s.png"%fnameTag)


