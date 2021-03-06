function out = subjFileInfo(subj,expt)

out = [];
out.dataPath                = ['~/Google Drive/Research/tACS/tACS_ER_task/data/' expt '/s' num2str(subj) '/' ];
out.infoStimFileName        = [];
out.stimulationFileName     = [];
out.eegEncodingFileName     = [];
out.eegRetrievalInfoFileName= [];
out.eegRetrievalFileName    = [];

switch expt
    case 'tacs_enc'
        switch  subj
            case 0
                warning([' although data was collected for this subject, markers did not '...
                    'work. '])
            case 1
                warning('no subject')
            case 2
                warning('self-paced retrieval task. changed from subject 4 on')
                out.infoStimFileName        = '20151121103036_s2.info';
                out.stimulationFileName     = '20151121103035_s2.stim';
                out.eegEncodingFileName     = '20151121103036_s2.easy';
                out.eegRetrievalInfoFileName= '20151121104805_s2.info';
                out.eegRetrievalFileName    = '20151121104805_s2.easy';
            case 3
                warning('self-paced retrieval task. changed from subject 4 on')
                out.infoStimFileName        = '20151121132947_s3.info';
                out.stimulationFileName     = '20151121132946_s3.stim';
                out.eegEncodingFileName     = '20151121132947_s3.easy';
                out.eegRetrievalInfoFileName= '20151121134627_s3.info';
                out.eegRetrievalFileName    = '20151121134627_s3.easy';
            case 4
                out.infoStimFileName        = '20151124132824_s4.info';
                out.stimulationFileName     = '20151124132824_s4.stim';
                out.eegEncodingFileName     = '20151124132824_s4.easy';
                out.eegRetrievalInfoFileName= '20151124134534_s4.info';
                out.eegRetrievalFileName    = '20151124134534_s4.easy';
            case 5
                out.infoStimFileName        = '20151124152044_s5.info';
                out.stimulationFileName     = '20151124152043_s5.stim';
                out.eegEncodingFileName     = '20151124152044_s5.easy';
                out.eegRetrievalInfoFileName= '20151124153808_s5.info';
                out.eegRetrievalFileName    = '20151124153808_s5.easy';
            case 6
                out.infoStimFileName        = '20151128103023_s6.info';
                out.stimulationFileName     = '20151128103023_s6.stim';
                out.eegEncodingFileName     = '20151128103023_s6.easy';
                out.eegRetrievalInfoFileName= '20151128104644_s6.info';
                out.eegRetrievalFileName    = '20151128104644_s6.easy';
            case 7
                out.infoStimFileName        = '20151202092812_s7.info';
                out.stimulationFileName     = '20151202092812_s7.stim';
                out.eegEncodingFileName     = '20151202092812_s7.easy';
                out.eegRetrievalInfoFileName= '20151202094551_s7.info';
                out.eegRetrievalFileName    = '20151202094551_s7.easy';
            case 8
                out.infoFileName            = '20151202112009_s8.info';
                out.stimulationFileName     = '20151202112009_s8.stim';
                out.eegEncodingFileName     = '20151202112009_s8.easy';
                out.eegRetrievalInfoFileName= '20151202113931_s8.info';
                out.eegRetrievalFileName    = '20151202113931_s8.easy';
            case 9
                out.infoStimFileName        = '20151202132023_s9.info';
                out.stimulationFileName     = '20151202132023_s9.stim';
                out.eegEncodingFileName     = '20151202132023_s9.easy';
                out.eegRetrievalInfoFileName= '20151202133758_s9.info';
                out.eegRetrievalFileName    = '20151202133758_s9.easy';
            case 10
                out.infoStimFileName        = '20151203142200_s10.info';
                out.stimulationFileName     = '20151203142200_s10.stim';
                out.eegEncodingFileName     = '20151203142200_s10.easy';
                out.eegRetrievalInfoFileName= '20151203144126_s10.info';
                out.eegRetrievalFileName    = '20151203144126_s10.easy';
            case 11
                out.infoStimFileName        = '20151206122419_s11.info';
                out.stimulationFileName     = '20151206122419_s11.stim';
                out.eegEncodingFileName     = '20151206122419_s11.easy';
                out.eegRetrievalInfoFileName= '20151206123950_s11.info';
                out.eegRetrievalFileName    = '20151206123950_s11.easy';
            otherwise
                error('subject does not exist')
        end
        
    case 'tacs_er_objstim'
        switch  subj
            case 1
                out.infoStimFileName        = '20160914114744_s1_enc.info';
                out.eegEncodingFileName     = '20160914114744_s1_enc.easy';
                out.stimulationFileName     = '20160914114456_s1_enc.stim';
                out.eegRetrievalInfoFileName= [];
                out.eegRetrievalFileName    = [];
            case 2
                out.infoStimFileName        = '20160914133746_s2_enc.info';
                out.eegEncodingFileName     = '20160914133746_s2_enc.easy';
                out.stimulationFileName     = '20160914133628_s2_enc.stim';
                out.eegRetrievalInfoFileName= '20160914140143_s2_ret.info';
                out.eegRetrievalFileName    = '20160914140143_s2_ret.easy';
            case 3
                out.infoStimFileName        = '20160920094433_s3_enc.info';
                out.stimulationFileName     = '20160920094210_s3_enc.stim';
                out.eegEncodingFileName     = '20160920094433_s3_enc.easy';
                out.eegRetrievalInfoFileName= '20160920101002_s3_ret.info';
                out.eegRetrievalFileName    = '20160920101002_s3_ret.easy';
            case 4
                out.infoStimFileName        = '20160920113551_s4_enc.info';
                out.stimulationFileName     = '20160920113344_s4_enc.stim';
                out.eegEncodingFileName     = '20160920113551_s4_enc.easy';
                out.eegRetrievalInfoFileName= '20160920120142_s4_ret.info';
                out.eegRetrievalFileName    = '20160920120142_s4_ret.easy';
            case 5
                out.infoStimFileName        = '20160920135635_s5_enc.info';
                out.stimulationFileName     = '20160920135513_s5_enc.stim';
                out.eegEncodingFileName     = '20160920135635_s5_enc.easy';
                out.eegRetrievalInfoFileName= '20160920142113_s5_ret.info';
                out.eegRetrievalFileName    = '20160920142113_s5_ret.easy';
            case 6
                out.infoStimFileName        = '20160922164204_s6_enc.info';
                out.stimulationFileName     = '20160922163920_s6_enc.stim';
                out.eegEncodingFileName     = '20160922164204_s6_enc.easy';
                out.eegRetrievalInfoFileName= '20160922170632_s6_ret.info';
                out.eegRetrievalFileName    = '20160922170632_s6_ret.easy';
            case 7
                out.infoStimFileName        = '20160929094117_s7_enc.info';
                out.stimulationFileName     = '20160929093943_s7_enc.stim';
                out.eegEncodingFileName     = '20160929094117_s7_enc.easy';
                out.eegRetrievalInfoFileName= '20160929100831_s7_ret.info';
                out.eegRetrievalFileName    = '20160929100831_s7_ret.easy';
            case 8
                out.infoStimFileName        = '20160929140157_s8_enc.info';
                out.stimulationFileName     = '20160929135923_s8_enc.stim';
                out.eegEncodingFileName     = '20160929140157_s8_enc.easy';
                out.eegRetrievalInfoFileName= '20160929142837_s8_ret.info';
                out.eegRetrievalFileName    = '20160929142837_s8_ret.easy';
            case 9
                out.infoStimFileName        = '20161003140732_s9_enc.info';
                out.stimulationFileName     = '20161003140600_s9_enc.stim';
                out.eegEncodingFileName     = '20161003140732_s9_enc.easy';
                out.eegRetrievalInfoFileName= '20161003143304_s9_ret.info';
                out.eegRetrievalFileName    = '20161003143304_s9_ret.easy';
            case 10
                out.infoStimFileName        = '20161003160702_s10_enc.info';
                out.stimulationFileName     = '20161003160523_s10_enc.stim';
                out.eegEncodingFileName     = '20161003160702_s10_enc.easy';
                out.eegRetrievalInfoFileName= '20161003163210_s10_ret.info';
                out.eegRetrievalFileName    = '20161003163210_s10_ret.easy';
            case 11
                out.infoStimFileName        = '20161005140259_s11_enc.info';
                out.stimulationFileName     = '20161005140151_s11_enc.stim';
                out.eegEncodingFileName     = '20161005140259_s11_enc.easy';
                out.eegRetrievalInfoFileName= '20161005142753_s11_ret.info';
                out.eegRetrievalFileName    = '20161005142753_s11_ret.easy';
            case 12
                out.infoStimFileName        = '20161006094742_s12_enc.info';
                out.stimulationFileName     = '20161006094500_s12_enc.stim';
                out.eegEncodingFileName     = '20161006094742_s12_enc.easy';
                out.eegRetrievalInfoFileName= '20161006101236_s12_ret.info';
                out.eegRetrievalFileName    = '20161006101236_s12_ret.easy';
            case 13
                out.infoStimFileName        = '20161006115009_s13_enc.info';
                out.stimulationFileName     = '20161006114854_s13_enc.stim';
                out.eegEncodingFileName     = '20161006115009_s13_enc.easy';
                out.eegRetrievalInfoFileName= '20161006121614_s13_ret.info';
                out.eegRetrievalFileName    = '20161006121614_s13_ret.easy';
            case 14
                out.infoStimFileName        = '20161006140737_s14_enc.info';
                out.stimulationFileName     = '20161006140429_s14_enc.stim';
                out.eegEncodingFileName     = '20161006140737_s14_enc.easy';
                out.eegRetrievalInfoFileName= '20161006143416_s14_ret.info';
                out.eegRetrievalFileName    = '20161006143416_s14_ret.easy';
            case 15
                out.infoStimFileName        = '20161010161737_s15_enc.info';
                out.stimulationFileName     = '20161010161534_s15_enc.stim';
                out.eegEncodingFileName     = '20161010161737_s15_enc.easy';
                out.eegRetrievalInfoFileName= '20161010164426_s15_ret.info';
                out.eegRetrievalFileName    = '20161010164426_s15_ret.easy';
            case 16
                out.infoStimFileName        = '20161011093521_s16_enc.info';
                out.stimulationFileName     = '20161011093352_s16_enc.stim';
                out.eegEncodingFileName     = '20161011093521_s16_enc.easy';
                out.eegRetrievalInfoFileName= '20161011100035_s16_ret.info';
                out.eegRetrievalFileName    = '20161011100035_s16_ret.easy';
            case 17
                out.infoStimFileName        = '20161011140200_s17_enc.info';
                out.stimulationFileName     = '20161011140024_s17_enc.stim';
                out.eegEncodingFileName     = '20161011140200_s17_enc.easy';
                out.eegRetrievalInfoFileName= '20161011142652_s17_ret.info';
                out.eegRetrievalFileName    = '20161011142652_s17_ret.easy';
            case 18
                out.infoStimFileName        = '20161011160500_s18_enc.info';
                out.stimulationFileName     = '20161011160338_s18_enc.stim';
                out.eegEncodingFileName     = '20161011160500_s18_enc.easy';
                out.eegRetrievalInfoFileName= '20161011163237_s18_ret.info';
                out.eegRetrievalFileName    = '20161011163237_s18_ret.easy';
            case 19
                out.infoStimFileName        = '20161019094331_s19_enc.info';
                out.stimulationFileName     = '20161019094206_s19_enc.stim';
                out.eegEncodingFileName     = '20161019094331_s19_enc.easy';
                out.eegRetrievalInfoFileName= '20161019100835_s19_ret.info';
                out.eegRetrievalFileName    = '20161019100835_s19_ret.easy';
            case 20
                out.infoStimFileName        = '20161019140313_s20_enc.info';
                out.stimulationFileName     = '20161019140059_s20_enc.stim';
                out.eegEncodingFileName     = '20161019140313_s20_enc.easy';
                out.eegRetrievalInfoFileName= '20161019142807_s20_ret.info';
                out.eegRetrievalFileName    = '20161019142807_s20_ret.easy';
            case 21
                out.infoStimFileName        = '20161023104355_s21_enc.info';
                out.stimulationFileName     = '20161023104154_s21_enc.stim';
                out.eegEncodingFileName     = '20161023104355_s21_enc.easy';
                out.eegRetrievalInfoFileName= '20161023110845_s21_ret.info';
                out.eegRetrievalFileName    = '20161023110845_s21_ret.easy';
            case 22
                out.infoStimFileName        = '20161023125554_s22_enc.info';
                out.stimulationFileName     = '20161023125335_s22_enc.stim';
                out.eegEncodingFileName     = '20161023125554_s22_enc.easy';
                out.eegRetrievalInfoFileName= '20161023132233_s22_ret.info';
                out.eegRetrievalFileName    = '20161023132233_s22_ret.easy';
            case 23
                out.infoStimFileName        = '20161027160939_s23_enc.info';
                out.stimulationFileName     = '20161027160626_s23_enc.stim';
                out.eegEncodingFileName     = '20161027160939_s23_enc.easy';
                out.eegRetrievalInfoFileName= '20161027163429_s23_ret.info';
                out.eegRetrievalFileName    = '20161027163429_s23_ret.easy';
            case 24
                out.infoStimFileName        = '20161029124508_s24_enc.info';
                out.stimulationFileName     = '20161029124313_s24_enc.stim';
                out.eegEncodingFileName     = '20161029124508_s24_enc.easy';
                out.eegRetrievalInfoFileName= '20161029131149_s24_ret.info';
                out.eegRetrievalFileName    = '20161029131149_s24_ret.easy';
            case 25
                out.infoStimFileName        = '20161029151715_s25_enc.info';
                out.stimulationFileName     = '20161029151448_s25_enc.stim';
                out.eegEncodingFileName     = '20161029151715_s25_enc.easy';
                out.eegRetrievalInfoFileName= '20161029154321_s25_ret.info';
                out.eegRetrievalFileName    = '20161029154321_s25_ret.easy';
            case 26
                out.infoStimFileName        = '20161120103658_s26_enc.info';
                out.stimulationFileName     = '20161120103451_s26_enc.stim';
                out.eegEncodingFileName     = '20161120103658_s26_enc.easy';
                out.eegRetrievalInfoFileName= '20161120110242_s26_ret.info';
                out.eegRetrievalFileName    = '20161120110242_s26_ret.easy';
            case 27
                out.infoStimFileName        = '20161120124713_s27_enc.info';
                out.stimulationFileName     = '20161120124501_s27_enc.stim';
                out.eegEncodingFileName     = '20161120124713_s27_enc.easy';
                out.eegRetrievalInfoFileName= '20161120131356_s27_ret.info';
                out.eegRetrievalFileName    = '20161120131356_s27_ret.easy';
            case 28
                out.infoStimFileName        = '20161121094552_s28_enc.info';
                out.stimulationFileName     = '20161121094300_s28_enc.stim';
                out.eegEncodingFileName     = '20161121094552_s28_enc.easy';
                out.eegRetrievalInfoFileName= '20161121101053_s28_ret.info';
                out.eegRetrievalFileName    = '20161121101053_s28_ret.easy';
            case 29
                out.infoStimFileName        = '20161121140358_s29_enc.info';
                out.stimulationFileName     = '20161121140148_s29_enc.stim';
                out.eegEncodingFileName     = '20161121140358_s29_enc.easy';
                out.eegRetrievalInfoFileName= '20161121142834_s29_ret.info';
                out.eegRetrievalFileName    = '20161121142834_s29_ret.easy';
            case 30
                out.infoStimFileName        = '20161121161230_s30_enc.info';
                out.stimulationFileName     = '20161121161012_s30_enc.stim';
                out.eegEncodingFileName     = '20161121161230_s30_enc.easy';
                out.eegRetrievalInfoFileName= '20161121163640_s30_ret.info';
                out.eegRetrievalFileName    = '20161121163640_s30_ret.easy';
            case 31
                out.infoStimFileName        = '20161121182114_s31_enc.info';
                out.stimulationFileName     = '20161121181906_s31_enc.stim';
                out.eegEncodingFileName     = '20161121182114_s31_enc.easy';
                out.eegRetrievalInfoFileName= '20161121184518_s31_ret.info';
                out.eegRetrievalFileName    = '20161121184518_s31_ret.easy';
            case 32
                out.infoStimFileName        = '20161122120046_s32_enc.info';
                out.stimulationFileName     = '20161122115750_s32_enc.stim';
                out.eegEncodingFileName     = '20161122120046_s32_enc.easy';
                out.eegRetrievalInfoFileName= '20161122122703_s32_ret.info';
                out.eegRetrievalFileName    = '20161122122703_s32_ret.easy';
            case 33
                out.infoStimFileName        = '20161129140416_s33_enc.info';
                out.stimulationFileName     = '20161129140200_s33_enc.stim';
                out.eegEncodingFileName     = '20161129140416_s33_enc.easy';
                out.eegRetrievalInfoFileName= '20161129142923_s33_ret.info';
                out.eegRetrievalFileName    = '20161129142923_s33_ret.easy';
            case 34
                out.infoStimFileName        = '20161129161533_s34_enc.info';
                out.stimulationFileName     = '20161129161332_s34_enc.stim';
                out.eegEncodingFileName     = '20161129161533_s34_enc.easy';
                out.eegRetrievalInfoFileName= '20161129164332_s34_ret.info';
                out.eegRetrievalFileName    = '20161129164332_s34_ret.easy';
            case 35
                out.infoStimFileName        = '20161130094051_s35_enc.info';
                out.stimulationFileName     = '20161130093833_s35_enc.stim';
                out.eegEncodingFileName     = '20161130094051_s35_enc.easy';
                out.eegRetrievalInfoFileName= '20161130100543_s35_ret.info';
                out.eegRetrievalFileName    = '20161130100543_s35_ret.easy';
            case 36
                out.infoStimFileName        = '20161130135323_s36_enc.info';
                out.stimulationFileName     = '20161130135039_s36_enc.stim';
                out.eegEncodingFileName     = '20161130135323_s36_enc.easy';
                out.eegRetrievalInfoFileName= '20161130142027_s36_ret.info';
                out.eegRetrievalFileName    = '20161130142027_s36_ret.easy';
            case 37
                out.infoStimFileName        = '20161130160858_s37_enc.info';
                out.stimulationFileName     = '20161130160613_s37_enc.stim';
                out.eegEncodingFileName     = '20161130160858_s37_enc.easy';
                out.eegRetrievalInfoFileName= '20161130163536_s37_ret.info';
                out.eegRetrievalFileName    = '20161130163536_s37_ret.easy';
        end
        
end

