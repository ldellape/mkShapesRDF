Steps = {
    "Summer23BPix_130x_nAODv12": {
        "isChain": True,
        "do4MC": True,
        "do4Data": True,
        "selection": '"((nElectron+nMuon)>0)"',
        "subTargets": [
            #"lumiMask",
            "leptonMaker",
            "fatjetSel",
            "FatJMECalculator",
            #"lepSel",
            #"jetSelUL",
            # "CleanFatJet",
            # "rochesterDATA",
            #"l2Kin",
            # "l3Kin",
            # "l4Kin",
            # "trigData",
            # "formulasDATA",
            "finalSnapshot_MC",
        ],
    },
    "Run2018_UL2018" : { 
        "isChain": True,
        "do4MC": True,
        "do4Data": True,
        "selection": '"((nElectron+nMuon)>0)"',
        "subTargets": [
            "lumiMask",
            "leptonMaker",
            "fatjetSel",
            #"lepSel",
            #"jetSelUL",
            # "CleanFatJet",
            # "rochesterDATA",
            #"l2Kin",
            # "l3Kin",
            # "l4Kin",
            # "trigData",
            # "formulasDATA",
            "finalSnapshot_DATA",
        ],
    },
    "JES_18": {
        "isChain": True,
        "do4MC": True,
        "do4Data": False,
        "subTargets": [
            "JES_modules_18UL",
            "l2Kin",
            "finalSnapshot_Variations",
        ],
    },
    "JES_18_test": {
        "isChain": True,
        "do4MC": True,
        "do4Data": False,
        "subTargets": ["JES_modules_18UL", "l2Kin", "histogram"],
    },
    "MCl1loose2018v9": {
        "isChain": True,
        "do4MC": True,
        "do4Data": False,
        "selection": '"((nElectron+nMuon)>0)"',
        "subTargets": [
           "leptonMaker",
           #"fatjetSel",
            "lepSel",
            "jetSelUL",
            "PromptParticlesGenVars",
            "GenVar",
            "GenLeptonMatch",
            "HiggsGenVars",
            "TopGenVars",
            "WGammaStar",
            "DressedLeptons",
        ],
        # 'wwNLL','ggHTheoryUncertainty', 'qqHTheoryUncertainty', 'EFTGen'],
    },
    "MCCorr2018v9": {
        "isChain": True,
        "do4MC": True,
        "do4Data": False,
        "subTargets": [
            "baseW",
            # "JES_modules_18UL",
            # "JERsMCUL",
            # # "FatJERsMCUL",
            "btagPerJet_DeepCSV_2018UL",
            "btagPerJet_DeepJet_2018UL",
            # "JetPUID_SF_UL",
            # "rochesterMC",
            # "trigMC",
            # "leptonSF",
            # "puW",
            "l2Kin",
            # # "l3Kin",
            # # "l4Kin",
            # "formulasMC",
            # # "EmbeddingVeto",
            # # "wwNLOEWK",
            # # "wwNLOEWK2",
            # # "wzNLOEWK",
            # # "zzNLOEWK",
            # # "zNLOEWK",
            # # "wNLOEWK",
            # # "qqHTheoryUncertainty",
            # # "CleanFatJet",
            # # "BoostedWtagSF",
            # "leptonMVAFiller",
        ],
    },
    "MCFull2018v9": {
        "isChain": True,
        "do4MC": True,
        "do4Data": False,
        "subTargets": [
            "MCl1loose2018v9",
            "MCCorr2018v9",
            "finalSnapshot_MC",
        ],
    },
    "leptonMaker": {
        "isChain": False,
        "do4MC": True,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.LeptonMaker",
        "declare": "leptonMaker = lambda : LeptonMaker()",
        "module": "leptonMaker()",
    },
    "lepSel": {
        "isChain": False,
        "do4MC": True,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.LeptonSel",
        "declare": 'leptonSel = lambda : LeptonSel("Loose", 1)',
        "module": "leptonSel()",
    },
    "jetSelUL": {
        "isChain": False,
        "do4MC": True,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.JetSel",
        # jetid=2,pujetid='loose',minpt=15.0,maxeta=4.7, UL2016fix=False "
        "declare": 'jetSel = lambda : JetSel(2,"loose",15.0,4.7,False)',
        "module": "jetSel()",
    },
    "lumiMask": {
        "isChain": False,
        "do4MC": False,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.LumiMask",
        "declare": "lumiMask = lambda : LumiMask(lumiFile)",
        "module": "lumiMask()",
    },
    "fatjetSel": {
        "isChain":False,
        "do4MC":True,
        "do4Data":True,
        "import": "mkShapesRDF.processor.modules.FatJetSel",
        "declare": "fatJetSel = lambda : FatJetSel()",
        "module" : "fatJetSel()",    
    },
    "FatJMECalculator": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.FatJMECalculatorRun3",
        #"declare": 'fatjmeCalculator = lambda : FatJMECalculator("RPLME_FW/processor/data/jsonpog-integration/POG/JME/2022_Summer22EE/fatJet_jerc.json.gz", "Summer22EE_22Sep2023_V2_MC", "Summer22EE_22Sep2023_JRV1_MC", "RPLME_FW/processor/data/jsonpog-integration/POG/JME/jer_smear.json.gz", "AK8PFPuppi", "RPLME_FW/processor/data/jsonpog-integration/POG/JME/2022_Summer22EE/jet_jerc.json.gz","Summer22EE_22Sep2023_V2_MC", "AK4PFPuppi",  do_JER=True, store_nominal=True, store_variations=True)',
        "declare": 'fatjmeCalculator = lambda : FatJMECalculator("RPLME_FW/processor/data/jsonpog-integration/POG/JME/2023_Summer23BPix/fatJet_jerc.json.gz", "Summer23BPixPrompt23_V1_MC", "Summer23BPixPrompt23_RunD_JRV1_MC", "RPLME_FW/processor/data/jsonpog-integration/POG/JME/jer_smear.json.gz", "AK8PFPuppi", "RPLME_FW/processor/data/jsonpog-integration/POG/JME/2023_Summer23BPix/jet_jerc.json.gz","Summer23BPixPrompt23_V1_MC", "AK4PFPuppi",  do_JER=True, store_nominal=True, store_variations=True)',
        "module": "fatjmeCalculator()",
    },
    "PromptParticlesGenVars": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.PromptParticlesGenVarsProducer",
        "declare": "PromptParticlesGenVars = lambda : PromptParticlesGenVarsProducer()",
        "module": "PromptParticlesGenVars()",
    },
    "GenVar": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.GenVarProducer",
        "declare": "GenVar = lambda : GenVarProducer()",
        "module": "GenVar()",
    },
    "GenLeptonMatch": {
        "isChain": False,
        "do4MC": True,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.GenLeptonMatchProducer",
        "declare": "GenLeptonMatch = lambda : GenLeptonMatchProducer()",
        "module": "GenLeptonMatch()",
    },
    "HiggsGenVars": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.HiggsGenVarsProducer",
        "declare": "HiggsGenVars = lambda : HiggsGenVarsProducer()",
        "module": "HiggsGenVars()",
    },
    "TopGenVars": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.TopGenVarsProducer",
        "declare": "TopGenVars = lambda : TopGenVarsProducer()",
        "module": "TopGenVars()",
    },
    "WGammaStar": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.WGammaStarProducer",
        "declare": "WGammaStar = lambda : WGammaStarProducer()",
        "module": "WGammaStar()",
    },
    "DressedLeptons": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.DressedLeptonProducer",
        "declare": "DressedLeptons = lambda : DressedLeptonProducer(0.3)",
        "module": "DressedLeptons()",
    },
    "baseW": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.BaseW",
        "declare": "baseW = lambda : BaseW(sampleName, files, xs_db, RPLME_genEventSumw)",
        "module": "baseW()",
    },
    "JES_modules_18UL": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.JMECalculator",
        "declare": 'jmeCalculator = lambda : JMECalculator("Summer19UL18_V5_MC", "Summer19UL18_JRV2_MC", \
            jet_object="AK4PFchs", do_Jets=True, do_MET=True, do_JER=False, store_nominal=False, store_variations=True)',
        "module": "jmeCalculator()",
    },
    "l2Kin": {
        "isChain": False,
        "do4MC": True,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.l2KinProducer",
        "declare": "l2Kin = lambda : l2KinProducer()",
        "module": "l2Kin()",
    },
    "leptonSF": {
        "isChain": False,
        "do4MC": False,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.LeptonSF",
        "declare": "leptonSF = lambda : LeptonSF(('RPLME_FW/processor/data/scale_factor/Full2018v9/electron.json.gz'))",
        "module": "leptonSF()",
    },
    "btagPerJet_DeepCSV_2018UL": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.btagSFProducerLatinos",
        "declare": 'btagPerJet_DeepCSV_2018UL = lambda : btagSFProducerLatinos(2018, "deepCSV", ["shape"], "shape", "RPLME_FW/processor/data/jsonpog-integration/POG/BTV/2018_UL/btagging.json.gz", ["jes","jesAbsolute","jesAbsolute_2018","jesBBEC1","jesBBEC1_2018","jesEC2","jesEC2_2018","jesFlavorQCD","jesHF","jesHF_2018","jesRelativeBal","jesRelativeSample_2018"])',
        "module": "btagPerJet_DeepCSV_2018UL()",
    },
    "btagPerJet_DeepJet_2018UL": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.btagSFProducerLatinos",
        "declare": 'btagPerJet_DeepJet_2018UL = lambda : btagSFProducerLatinos(2018, "deepJet", ["shape"], "shape", "RPLME_FW/processor/data/jsonpog-integration/POG/BTV/2018_UL/btagging.json.gz", ["jes","jesAbsolute","jesAbsolute_2018","jesBBEC1","jesBBEC1_2018","jesEC2","jesEC2_2018","jesFlavorQCD","jesHF","jesHF_2018","jesRelativeBal","jesRelativeSample_2018"])',
        "module": "btagPerJet_DeepJet_2018UL()",
    },
    "finalSnapshot_MC": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.Snapshot",
        "declare": "snapshot = lambda : Snapshot( \
                tmpOutputFilename='output.root', \
                columns=['*'], \
                eosPath='RPLME_EOSPATH', outputFilename='RPLME_OUTPUTFILENAME', \
                includeVariations=False, splitVariations=False, storeNominals=True )",
        "module": "snapshot()",
    },
    "finalSnapshot_Variations": {
        "isChain": False,
        "do4MC": True,
        "do4Data": False,
        "import": "mkShapesRDF.processor.modules.Snapshot",
        "declare": "snapshot = lambda : Snapshot( \
                tmpOutputFilename='output.root', \
                columns=['*'], \
                eosPath='RPLME_EOSPATH', outputFilename='RPLME_OUTPUTFILENAME', \
                includeVariations=True, splitVariations=True, storeNominals=False )",
        "module": "snapshot()",
    },
    "finalSnapshot_DATA": {
        "isChain": False,
        "do4MC": False,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.Snapshot",
        "declare": "snapshot = lambda : Snapshot( \
                tmpOutputFilename='output.root', \
                columns=['*'], \
                eosPath='RPLME_EOSPATH', outputFilename='RPLME_OUTPUTFILENAME', \
                includeVariations=False, splitVariations=False, storeNominals=True )",
        "module": "snapshot()",
    },
    "histogram": {
        "isChain": False,
        "do4MC": True,
        "do4Data": True,
        "import": "mkShapesRDF.processor.modules.Histogram",
        "declare": "histogram = lambda : Histogram( \
                outputFilename='output2.root', \
                variables=[(('mjj', 'mjj', 100, 10, 1000), 'new_fw_mjj', 'baseW * genWeight')] \
                    )",
        "module": "histogram()",
    },
}
