if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    from CRABClient.UserUtilities import config, getUsernameFromSiteDB
    config = config()

    config.General.workArea = 'VertexCompositeAna'
    config.General.transferOutputs = True
    config.General.transferLogs = False
    config.JobType.pluginName = 'Analysis'
    config.JobType.maxMemoryMB = 3000
    config.JobType.maxJobRuntimeMin = 2750
#    config.JobType.psetName = '../test/pPbFlowCorrSkim_2016_D0_cfg.py'
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 1000
    config.Data.splitting = 'FileBased'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
    config.Data.publication = True
    config.Data.inputDBS = 'phys03'
    config.Site.storageSite = 'T2_US_MIT'
#    config.Site.storageSite = 'T3_US_Rice'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'PbPb2018SkimPv1'
    config.JobType.psetName = '../test/PbPbSkim2018_D0Both_mc_cfg.py'
    config.Data.outputDatasetTag = 'prompt_pt1p2_y2p4_hi1031p1_Skim_v1'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-prompt_pt1p2_y2p4_hi1031p1_RECO_v1-4ec1ff910eab1786049c0e4c444e87c6/USER'
    submit(config)

    config.General.requestName = 'PbPb2018SkimNPv1'
    config.Data.outputDatasetTag = 'nonprompt_pt1p2_y2p4_hi1031p1_Skim_v2'
    config.Data.inputDataset = '/Pythia8_TuneCUETP8M1_13TeV_D0_PiK/davidlw-nonprompt_pt1p2_y2p4_hi1031p1_RECO_v2-4ec1ff910eab1786049c0e4c444e87c6/USER'
    submit(config)
