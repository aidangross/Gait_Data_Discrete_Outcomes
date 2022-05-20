close all
clear 
clc


test = process('S1');
function y = process(sheet)
prompt = 'Enter Study ID';
studyid = input(prompt,'s');
prompt = 'Enter Shoe Condition';
shoe_type = input(prompt,'s');
[file_list,path_n] = uigetfile('.txt','Select files for processing','MultiSelect','on');
if iscell(file_list) == 0
    file_list = {file_list};
end

y_total = zeros(1,88);

for i = 1:length(file_list)
    filename = file_list{i};
    data = readmatrix([path_n filename]);
    height = data(1,2);
    mass = data(1,3);
    speed = data(1,4);
    
    %Ground reaction force plate 2:
    fp2 = data(:,6);
    fp2 = (rmmissing(fp2))/9.81;
    fifty = round(0.5*round(length(fp2)));
    
    %Peak ground reaction forces for first and second peaks:
    %1
    [grf1, grf1_ind] = max(fp2(1:fifty));
    %2
    [grf2, grf2_ind] = max(fp2(fifty:end));
    
    %Normalizing GRF to %BW
    grf1 = (grf1/mass)*100;
    grf2 = (grf2/mass)*100; 
    
    %Medial Lateral GRF:
    medial_lateral_grf = data(:,62);
    medial_lateral_grf = (rmmissing(medial_lateral_grf))/9.81;
    
    %Peak Medial GRF:
    p_m_grf = (min(medial_lateral_grf/mass))*100;
    
    %Peak Lateral GRF: 
    p_l_grf = (max(medial_lateral_grf/mass))*100;
    
    %%%%%%%%%%%%%%%%%
    %%%%% KNEE %%%%%%
    %%%%%%%%%%%%%%%%%
   
    %Knee flexion angle in first frame of stance - heel strike:  
    kneex = data(:,55);
    kneex = rmmissing(kneex);
    %Knee flexion at heel strike:
    kfhs = kneex(1);
    %Knee flexion at toe off:
    kfto = kneex(end);
   
    %Local maximum of knee flexion angle in first 40% of stance:
    mskf = min(kneex(1:round(.4*length(kneex))));
    
    %Local minimum of knee flexion in second half of stance:
    mske_peak_stance = max(kneex(round(.5*length(kneex)):end));
    
    %Knee flexion range of motion for stance phase:
    kfrom_stance = max(kneex)-min(kneex);
    
    %Knee rotation angle:
    kneez = data(:,57);
    kneez = rmmissing(kneez);
    %Knee rotation angle:
    krot_avg = mean(kneez);
    %Knee rotation at heel strike:
    krot_hs = kneez(1);
    
    %Knee adduction/abduction angle:
    kneey = data(:,56);
    kneey = rmmissing(kneey);
    [peak_varus,p_stance_peak_varus] = max(kneey);
   
    %%%%%%%%%%%%%%%%%
    %%%%% HIP %%%%%%%
    %%%%%%%%%%%%%%%%%
    
    %Hip flexion angle in stance phase:
    hipx_stance = data(:,58);
    hipx_stance = rmmissing(hipx_stance);
    %Hip flexion at heel strike:
    hfhs = hipx_stance(1);
    %Hip flexion toe off:
    hfto = hipx_stance(end);
    
    %Max hip flexion in first 20% of stance:
    hf_max = max(hipx_stance(1:round(.2*length(hipx_stance))));
    %Max hip extension = local minimum in final 30% of stance:
    he_max = min(hipx_stance(round(.7*length(hipx_stance)):end));
    
    %Range of motion of hip angle during stance:
    hsrom_stance = hf_max-he_max;
    
    %%%%%%%%%%%%%%%%%
    %%%%% ANKLE %%%%%
    %%%%%%%%%%%%%%%%%
    
    %Ankle flexion angle stance phase:
    anklex_stance = data(:,40);
    anklex_stance = rmmissing(anklex_stance);
    %Ankle flexion angle at heel strike:
    afhs = anklex_stance(1);
    %Ankle flexion angle at toe off:
    afto = anklex_stance(end);
    %Minimum angle in stance (assuming pf is negative):
    apf_peak = min(anklex_stance);
    %Max positive value - max negative value in stance: 
    asrom_stance = max(anklex_stance)-min(anklex_stance);
    
    %Frontal plane ankle angle stance:
    ankley = data(:,41);
    ankley = rmmissing(ankley);
   
    %Ankle inversion/eversion:
    [aev_peak,p_aev_peak] = min(ankley);
    [aiv_peak,p_aiv_peak] = max(ankley);
    mean_ankle_y = mean(ankley);
    
    %Ankle inversion/eversion peak velocity in early stance: 
    aie_p_v = max(abs(diff(ankley(1:round(.5*length(ankley))))));
    
    %%%%%%%%%%%%%%%%%%%
    %%%%% PELVIS %%%%%%
    %%%%%%%%%%%%%%%%%%%
    
    %Anterior - Posterior pelvis movement stance phase (sagittal plane):
    pelvisxs = data(:,7);
    pelvisxs = rmmissing(pelvisxs);
    pant_peak = min(pelvisxs);
    ppost_peak = max (pelvisxs);
    
    %Pelvis Hike - Drop movement (frontal plane):
    pelvisy = data(:,8);
    pelvisy = rmmissing(pelvisy);
    pdrop_peak = min(pelvisy);  
    phike_peak = max(pelvisy);
    
    %Pelvis forwards/backwards rotation (transverse):
    pelvisz = data(:,9);
    pelvisz = rmmissing(pelvisz);
    prot_peak = max(pelvisz);
    
    %%%%%%%%%%%%%%%%%%
    %%%KNEE MOMENTS%%%
    %%%%%%%%%%%%%%%%%%
    
    %Knee flexion moment:
    kfm = data(:,13);
    kfm = rmmissing(kfm);
 
    %Knee flexion/extension peak moments: 
    first_kem_peak = min(kfm(1:round(.2*length(kfm))));
    kfm_peak = max(kfm(1:round(.5*length(kfm))));
    second_kem_peak = min(kfm(round(.5*length(kfm)):end));
    
    %Knee adduction moment:
    kam = data(:,14);
    kam = rmmissing(kam);
   
    %Knee adduction moment peaks:
    kam_first = min(kam(1:round(.5*length(kam))));
    kam_second = min(kam(round(.5*length(kam)):end));
    
    %Knee rotation moment:
    krot = data(:,15);
    krot = rmmissing(krot);
    %Minimum knee rotation moment:
    kirot = min(krot);
    %Maximum knee rotation moment:
    kerot = max(krot);
    
    %%%%%%%%%%%%%%%%%
    %%%HIP MOMENTS%%%
    %%%%%%%%%%%%%%%%%
    
    %Hip flexion moment:
    hipm = data(:,10);
    hipm = rmmissing(hipm);
   
    %Max hip flexiom moment:
    hfm_peak = min(hipm(1:round(.5*length(hipm))));
    %Max hip extension moment:
    hem_peak = max(hipm(round(.5*length(hipm)):end));
    
    %%%%%%%%%%%%%%%%%
    %%Ankle MOMENTS%%
    %%%%%%%%%%%%%%%%%
    
    %Ankle flexion moment: 
    ank_flex_mom = data(:,16);
    ank_flex_mom = rmmissing(ank_flex_mom);
    %Local max first 50% stance
    adm_peak = max(ank_flex_mom(1:round(.5*length(ank_flex_mom))));
    %Local min second 50% stance
    apm_peak = min(ank_flex_mom(round(.5*length(ank_flex_mom)):end));
    
    %Ankle inversion/eversion moment:
    aiem = data(:,17);
    aiem = rmmissing(aiem);
   
    %ankle inversion moment peaks:
    aim_first = min(aiem(1:round(.5*length(aiem))));
    aim_second = min(aiem(round(.5*length(aiem)):end));
    
    %ankle eversion moment peak:
    aem_peak = max(aiem(round(.25*length(aiem)):(.75*length(aiem))));
   
    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%Power calculations%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%
    hip_power= data(:,19); 
    knee_power= data(:,20);
    ankle_power= data(:,21);
    hip_sag_moment = data(:,10);
    peak_ext_pos_power = max(hip_power(hip_sag_moment < 0));
    peak_flex_pos_power = max(hip_power(hip_sag_moment > 0));
    hip_pos_power = max(hip_power);
    hip_neg_power = min(hip_power);
    knee_pos_power= max(knee_power);
    knee_neg_power= min(knee_power);
    ankle_pos_power= max(ankle_power);
    ankle_neg_power= min(ankle_power);
    total_pos_power= hip_pos_power+knee_pos_power+ankle_pos_power;
    total_neg_power= hip_neg_power+knee_neg_power+ankle_neg_power;
    hip_rel_pos_power = (hip_pos_power / total_pos_power) * 100;
    knee_rel_pos_power = (knee_pos_power / total_pos_power) * 100;
    ankle_rel_pos_power = (ankle_pos_power / total_pos_power) * 100;
    hip_rel_neg_power = (hip_neg_power / total_neg_power) * 100;
    knee_rel_neg_power = (knee_neg_power / total_neg_power) * 100;
    ankle_rel_neg_power = (ankle_neg_power / total_neg_power) * 100;
    frame_rate = 200;
    hip_pos_work = trapz(hip_power(hip_power > 0)) ./ frame_rate;
    hip_neg_work = trapz(hip_power(hip_power < 0)) ./ frame_rate;
    knee_pos_work = trapz(knee_power(knee_power > 0)) ./ frame_rate;
    knee_neg_work = trapz(knee_power(knee_power < 0)) ./ frame_rate;
    ankle_pos_work = trapz(ankle_power(ankle_power > 0)) ./ frame_rate;
    ankle_neg_work = trapz(ankle_power(ankle_power < 0)) ./ frame_rate;
    total_pos_work = hip_pos_work + knee_pos_work + ankle_pos_work;
    total_neg_work = hip_neg_work + knee_neg_work + ankle_neg_work;
    hip_rel_pos_work = (hip_pos_work / total_pos_work) * 100;
    knee_rel_pos_work = (knee_pos_work / total_pos_work) * 100;
    ankle_rel_pos_work = (ankle_pos_work / total_pos_work) * 100;
    hip_rel_neg_work = (hip_neg_work / total_neg_work) * 100;
    knee_rel_neg_work = (knee_neg_work / total_neg_work) * 100;
    ankle_rel_neg_work = (ankle_neg_work / total_neg_work) * 100;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Center of Pressure %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Medial Lateral GRF:
    medial_lateral_cop = data(:,61);
    medial_lateral_cop = rmmissing(medial_lateral_cop);
    
    %Peak Medial GRF:
    [p_m_cop, p_m_copi] = (max(medial_lateral_cop));
    
    %Peak Lateral GRF: 
    [p_l_cop, p_l_copi] = (min(medial_lateral_cop));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Loading Rate %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Calculate linear (LLR) and instantaneous loading rate (ILR)
    % ILR is integral of GRF divided by frame rate
     frame_rate = 2000;
     load_integral= trapz(fp2(1:fifty)) ./ frame_rate;
     ILR= max(load_integral);

     % LLR is peak GRF divided by time to peak
     [max_num,max_idx]= max(fp2(1:fifty));
     peak_time= max_idx ./ frame_rate;
     LLR= max_num ./ peak_time;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Multi-Segment Foot%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Rear-foot to Fore-foot: 
    cal_mety = data(:,29);
    cal_mety = rmmissing(cal_mety);
    
    %Mean calacaneus metatarsal angle over stance: 
    mean_cal_mety = mean(cal_mety);
    cal_mety_p_v = max(abs(diff(cal_mety(1:round(.5*length(cal_mety))))));
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     y = {studyid shoe_type mass height speed grf1 grf1_ind grf2 grf2_ind... 
      kfhs mskf mske_peak_stance kfrom_stance kfto mean_ankle_y...
      hfhs hfto hf_max he_max hsrom_stance afhs afto apf_peak asrom_stance...
      krot_avg krot_hs aev_peak p_aev_peak pant_peak pdrop_peak prot_peak first_kem_peak kfm_peak...
      second_kem_peak kam_first kam_second kerot kirot hfm_peak hem_peak apm_peak adm_peak... 
      peak_ext_pos_power peak_flex_pos_power hip_pos_power hip_neg_power...
      knee_pos_power knee_neg_power ankle_pos_power ankle_neg_power total_pos_power total_neg_power...
      hip_rel_pos_power knee_rel_pos_power ankle_rel_pos_power hip_rel_neg_power knee_rel_neg_power...
      ankle_rel_neg_power hip_pos_work hip_neg_work knee_pos_work knee_neg_work ankle_pos_work...
      ankle_neg_work total_pos_work total_neg_work hip_rel_pos_work knee_rel_pos_work...
      ankle_rel_pos_work hip_rel_neg_work knee_rel_neg_work ankle_rel_neg_work peak_varus...
      p_stance_peak_varus aim_first aim_second aiv_peak p_aiv_peak ILR LLR phike_peak...
      ppost_peak aem_peak p_m_grf p_l_grf p_m_cop p_m_copi p_l_cop p_l_copi aie_p_v mean_cal_mety...
      cal_mety_p_v};
     
     file = '/Users/aidangross/Desktop/Ind_Study/oofos_2_processing/OP2_Gait_Discrete_Variables.xlsx';
     
    
    %Conditional statements for correct sheets in excel document% 
    
    if contains(filename,'OoOG') 
    
        writecell(y,file,'Sheet','OP2_OoOG','WriteMode','append')
        
    elseif contains(filename,'OOmg')

        writecell(y,file,'Sheet','OP2_OOmg','WriteMode','append')
        
    elseif contains(filename,'Adidas')

        writecell(y,file,'Sheet','OP2_Adidas','WriteMode','append')
    
    elseif contains(filename,'Reef')

        writecell(y,file,'Sheet','OP2_Reef','WriteMode','append') 
    
    elseif contains(filename,'Birk')

        writecell(y,file,'Sheet','OP2_Birk','WriteMode','append') 
   
    elseif contains(filename,'Eva')

        writecell(y,file,'Sheet','OP2_Eva','WriteMode','append') 
    end 
  
    
    %Calculates the total for y over x number of trials
    y_total = y_total + cell2mat(y(1,5:92));
      
end

%Calculate the mean using the counter in the for loop(i.e.length(file_list))%

y_mean = y_total ./ length(file_list);

y_mean_cell = num2cell(y_mean);

y_mean_cell_labs = {studyid shoe_type mass height};

y_mean_cell_final = [y_mean_cell_labs y_mean_cell];

if contains(filename,'OoOG')    
    
    writecell(y_mean_cell_final,file,'Sheet','OP2_OoOG_Mean','WriteMode','append')
        
elseif contains(filename,'OOmg')

    writecell(y_mean_cell_final,file,'Sheet','OP2_OOmg_Mean','WriteMode','append')
    
elseif contains(filename,'Adidas')

    writecell(y_mean_cell_final,file,'Sheet','OP2_Adidas_Mean','WriteMode','append')

elseif contains(filename,'Reef')

    writecell(y_mean_cell_final,file,'Sheet','OP2_Reef_Mean','WriteMode','append')

elseif contains(filename,'Birk')

    writecell(y_mean_cell_final,file,'Sheet','OP2_Birk_Mean','WriteMode','append')
    
elseif contains(filename,'Eva')

    writecell(y_mean_cell_final,file,'Sheet','OP2_Eva_Mean','WriteMode','append')
end
end




