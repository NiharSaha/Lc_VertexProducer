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
#    config.JobType.maxMemoryMB = 3000
#    config.JobType.maxJobRuntimeMin = 2750
    config.JobType.psetName = '../test/pPbSkim2016_LambdaC_mc_cfg.py'
    config.Data.unitsPerJob = 1
#    config.Data.totalUnits = 1000
    config.Data.splitting = 'FileBased'
    config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#    config.Data.inputDBS = 'phys03'
    config.Data.publication = True
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

    config.General.requestName = 'pPb2016_pPbMC_Skim_LambdaCToKsP_v3'
    config.Data.outputDatasetTag = 'pPb_Skim_LambdaCToKsP_v3'
    config.Data.inputDataset = '/LambdaC-KsPr_LCpT-0p9_pPb-EmbEPOS_8p16_Pythia8/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v1/AODSIM'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMC_Skim_LambdaCToKsP_Pt6_v3'
    config.Data.outputDatasetTag = 'pPb_Skim_LambdaCToKsP_v3'
    config.Data.inputDataset = '/LambdaC-KsPr_LCpT-5p9_pPb-EmbEPOS_8p16_Pythia8/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v1/AODSIM'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMC_Skim_LambdaCToLamPi_v3'
    config.Data.outputDatasetTag = 'pPb_Skim_LambdaCToLamPi_v3'
    config.Data.inputDataset = '/LambdaC-LamPi_LCpT-0p9_pPb-EmbEPOS_8p16_Pythia8/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v1/AODSIM'
    submit(config)

    config.General.requestName = 'pPb2016_pPbMC_Skim_LambdaCToLamPi_Pt6_v3'
    config.Data.outputDatasetTag = 'pPb_Skim_LambdaCToLamPi_v3'
    config.Data.inputDataset = '/LambdaC-LamPi_LCpT-5p9_pPb-EmbEPOS_8p16_Pythia8/pPb816Summer16DR-pPbEmb_80X_mcRun2_pA_v4-v1/AODSIM'
    submit(config)
