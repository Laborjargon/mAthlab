function extractorMain_fm07(sub,const,scr,eyeana)
% ----------------------------------------------------------------------
% extractorMain_fm07(sub,const,scr,eyeana)
% ----------------------------------------------------------------------
% Goal of the function :
% Extract results, mean, number of the different combination of variable 
% used in this experiment.
% ----------------------------------------------------------------------
% Input(s) :
% sub   : subject configuration
% const : constant settings
% scr   : screen settings
% eyeana: eye analysis settings
% ----------------------------------------------------------------------
% Output(s):
% none
% ----------------------------------------------------------------------
% Function created by Nina HANNING (hanning.nina@gmail.com)
% Last update : 2025-01-14
% Project : noBlank
% Version : 2.0
% ----------------------------------------------------------------------

plot_0 = 1; % plot figure with timing info
plot_1 = 1; % plot perceptual reports


fileRes = csvread(sprintf('csv/%s_AllB_corMat.csv',sub.ini));
[~,col] = size(fileRes);



% load visVec 
% (t x nbf matrix containing information about target visibility (0 = blank, 1 = fully visible)
load(sprintf('csv/%s_AllB_visVecCor.mat',sub.ini));                  
visVec = visVecCor_AllB;


load(sprintf('trace/%s_AllB_corTrace.mat',sub.ini)); 
traceMat = corTrace_AllB;

if ~isdir('../../../Figures'); mkdir('../../../Figures'); end


%% Treat data
% =========================================================================


% trial info
block_col       =   1;      % block nb
trial_col       =   2;      % trial nb
cond_col        =   3;      % condition (1 = bg pulse, 2 = target pulse, 3 = stripes, 4 = classic)
st_col          =   7;      % saccade direction (left -1, right 1)
ecc_col         =   8;      % saccade direction (left -1, right 1)

gap_col         =   9;      % gap duration (0 or 250ms)
displ_col       =  10;      % target displacement (-2:2dva, values in pix, negative means leftward; divide by const.ppd to convert to dva)
pulse_jit_col   =  11;
stripe_jit_col  =  12;
stripe_dir_col  =  13;

blank_tt_col    =  14;      % planned blank type (1 = early jump /no blank, 2 = late jump /blank)

resp_col        =  21;      % response given (-1 left,1 right)
res_col         =  23;      % response correct (1,0)

% saccade parameters
trialON_col     =  51;      % tial sart time (el timestamp)
sacON_col       =  62;      % saccade onset time (el timestamp)
sacOFF_col      =  63;      % saccade offset time (el timestamp)
sacDur_col      =  64;      % saccade offset time (el timestamp)

sacAmp_col      =  67;      % saccade amplitude
sacxOnset       =  68;      % saccade onset x coordinate
sacyOnset       =  69;      % saccade onset y coordinate
sacxOffset      =  70;      % saccade offset x coordinate
sacyOffset      =  71;      % saccade offset y coordinate
sacLat_col      =  72;      % saccade latency
sacErr_col      =  73;      % saccade error (distance between offset and target center in dva)

tarPosX         =  74;      % target x coordinate
tarPosY         =  75;      % target y coordinate

resp2_col       = col+1;    % response given (1 right, 0 not right / left; like resp_col, just easier to average)
displ2_col      = col+2;    % response given (1 right, 0 not right / left; like resp_col, just easier to average)

sac_relON_col   = col+3;    % sacON relative to trial start (nbf)
sac_relOFF_col  = col+4;    % sacOFF relative to trial start (nbf)

sacONdurB_col   = col+5;    % sacON during blank (1 = yes, 2 = no)
sacOFFdurB_col  = col+6;    % sacOFF during blank (1 = yes, 2 = no)

blankType_col   = col+7;    % sacON and/or sacOFF during blank?
                            % 1 = sacON no, sacOFF no
                            % 2 = sacON no, sacOFF yes
                            % 3 = sacON yes, sacOFF yes
                            % 4 = sacON yes, sacOFF no

sacON_reBl_col   = col+ 8;  % sacON  relative to next blank
sacOFF_reBl_col  = col+ 9;  % sacOFF relative to next blank
sacON_reVis_col  = col+10;  % sacON  relative to next vis
sacOFF_reVis_col = col+11;  % sacOFF relative to next vis


shiftVal_col     = col+12;  % x-shift visVec to align trials for fig0
               


%% New computations

% separate condition 4 (classic blanking control)
% cond4idx = fileRes(:,cond_col)==4; % locate condition 4 trials
% 
% fileRes4 = fileRes(cond4idx,:); % matrix with only cond4 trials
% visVec4 = visVec(cond4idx,:);
% 
% fileRes = fileRes(~cond4idx,:); % matrix with all but cond4 trials
% visVec = visVec(~cond4idx,:);




%figure; % verify conversion
%subplot(3,1,1); plot(unique(visVec(~isnan(visVec)))');
%subplot(3,1,2); plot(visVec(1,:));
visVec(visVec>0 & visVec<1) = 0.5; % convert gradual values to 0.5
%subplot(3,1,3); plot(visVec(1,:));



% resp2_col, displ2_col
% convert -l/+r to 0inward/1outward
fileRes(:,resp2_col)  = fileRes(:,resp_col).*fileRes(:,st_col);     % flip response and displacement of leftward trials
fileRes(:,displ2_col) = fileRes(:,displ_col).*fileRes(:,st_col);
fileRes(fileRes(:,resp2_col)==-1,resp2_col) = 0;

fileRes(:,displ_col) = fileRes(:,displ_col)/const.ppd;  % convert displacement (pix -> dva)
fileRes(:,displ_col) = round(fileRes(:,displ_col),2);   


% sac_relON_col, sac_relOFF_col (frame nb relative to trial start)
fileRes(:,sac_relON_col)  = ((fileRes(:,sacON_col)  - fileRes(:,trialON_col))/1000 / scr.fd) - 1;
fileRes(:,sac_relOFF_col) = ((fileRes(:,sacOFF_col) - fileRes(:,trialON_col))/1000 / scr.fd) - 1;




% sacONdurB_col, sacOFFdurB_col
% blankType_col
% sacON_reBl_col, sacOFF_reBl_col
% sacON_reVis_col, sacOFF_reVis_col
fileRes(:,sacONdurB_col)    = NaN;
fileRes(:,sacOFFdurB_col)   = NaN;

fileRes(:,blankType_col)    = NaN;

fileRes(:,sacON_reBl_col)   = NaN;
fileRes(:,sacOFF_reBl_col)  = NaN;
fileRes(:,sacON_reVis_col)  = NaN;
fileRes(:,sacOFF_reVis_col) = NaN;
fileRes(:,shiftVal_col)     = NaN;

for i = 1:size(fileRes,1)
   
    nbf_sacON  = floor(fileRes(i,sac_relON_col));
    nbf_sacOFF = floor(fileRes(i,sac_relOFF_col));
        
    % sacON before blank (0) or during (1) ?
    if visVec(i,nbf_sacON)
        fileRes(i,sacONdurB_col) = 0;
    else
        fileRes(i,sacONdurB_col) = 1;
    end
    
    % sacOFF before blank (0) or during (1) ?
    if visVec(i,nbf_sacOFF)
        fileRes(i,sacOFFdurB_col) = 0;
    else
        fileRes(i,sacOFFdurB_col) = 1;
    end
    
    
    visVec_ON = visVec(i,nbf_sacON:end);
    visVec_OFF = visVec(i,nbf_sacOFF:end);
    
    
    %%% [ fm07-1 ]
    if fileRes(i,cond_col) == 4
        
        
        % condition 4 trials without blank should be categorized as 
        % 'blankType' 5
        % condition 4 trials with blank should be categorized as 
        % 'blankType' 6 
        
        % idea ben:
        if fileRes(i, gap_col) == 0
          fileRes(i,blankType_col) = 5;
        else
          fileRes(i,blankType_col) = 6;
        end

        % fileRes(i,blankType_col) = 
        
    else
        
        % [1] sacON no & sacOFF no
        % = classic 'no blank'
        if fileRes(i,sacONdurB_col) == 0 && fileRes(i,sacOFFdurB_col) == 0
            fileRes(i,blankType_col) = 1;
            
            fileRes(i,sacON_reBl_col)  = -find(visVec_ON==0,1);     % sacON  relative to next blank
            fileRes(i,sacOFF_reBl_col) = -find(visVec_OFF==0,1);    % sacOFF relative to next blank
            
        % [2] sacON no, sacOFF *BLANK*
        % = classic 'blanking'
        elseif fileRes(i,sacONdurB_col) == 0 && fileRes(i,sacOFFdurB_col) == 1
            fileRes(i,blankType_col) = 2;
            
            fileRes(i,sacON_reBl_col)   = -find(visVec_ON==0,1);    % sacON  relative to next blank
            fileRes(i,sacOFF_reVis_col) = -find(visVec_OFF~=0,1);   % sacOFF relative to next vis
            
        % [3] sacON *BLANK* & sacOFF *BLANK*
        % = (saccade during blank)
        elseif fileRes(i,sacONdurB_col) == 1 && fileRes(i,sacOFFdurB_col) == 1
            fileRes(i,blankType_col) = 3;
            
            fileRes(i,sacON_reVis_col)  = -find(visVec_ON~=0,1);    % sacON  relative to next vis
            fileRes(i,sacOFF_reVis_col) = -find(visVec_OFF~=0,1);   % sacOFF relative to next vis
            
        % [4] sacON *BLANK* & sacOFF no
        % = (blanking before saccade)
        elseif fileRes(i,sacONdurB_col) == 1 && fileRes(i,sacOFFdurB_col) == 0
            fileRes(i,blankType_col) = 4;
            
            fileRes(i,sacON_reVis_col) = -find(visVec_ON~=0,1);     % sacON  relative to next vis
            fileRes(i,sacOFF_reBl_col) = -find(visVec_OFF==0,1);    % sacOFF relative to next blank
        else
            error; % should not exist
        end
        
    end
    
end





%figure; hold on; % verify alignment
for i = 1:size(visVec,1)
    if fileRes(i,cond_col) ~= 4
        if fileRes(i,cond_col) == 3
            fileRes(i,shiftVal_col) = round(fileRes(i,stripe_jit_col)/const.stripeSp_ppf);
            %plot([1:size(visVec,2)]+fileRes(i,shiftVal_col),visVec(i,:));
        else
            
            fileRes(i,shiftVal_col) = -fileRes(i,pulse_jit_col);
            %plot([1:size(visVec,2)]+fileRes(i,shiftVal_col),visVec(i,:));
        end
    end
end


if plot_0
    
    fileRes_all = fileRes;
    
    fig0 = figure;
    
    colFade   = [170, 165, 57]./255;
    colBl     = [0,0,0];
    
    colEcc{1} = [36, 96, 104]./255;
    colEcc{2} = [69, 47, 116]./255;
    colEcc{3} = [147, 49, 87]./255;
   
    
    
    %---------------------------
    % [A] sacOn & sacOFF re blanking
    %---------------------------
    
    subplot(3,2,1:2); hold on;

    if const.condition == 3 % stripe condition
        title(sprintf('stripeW: %1.1f dva    |    stripeD: %1.1f dva    |    speed: %1.1f dps    |    fixRad: %1.3f dva',const.stripeW_deg, const.stripeD_deg, const.stripeSp_dps, const.fixRadVal));

        text(10,0.9,sprintf('durRamp:  %3.0f ms',const.durRamp*1000));
        text(10,0.8,sprintf('durBlank: %3.0f ms',const.durBlank*1000));
        text(10,0.7,sprintf('durSee:   %3.0f ms',const.durSee*1000));
        
    else % color pulse conditions
        title(sprintf('durBlank: %i ms    |    durSee: %i ms    |    durRamp: %i ms',const.durBlank*1000, const.durSee*1000, const.durRamp*1000));
    end
        
    % display stripes / pulses
    % match times
    
    % take 'longest' trial (and shift 'shorter' trials rightwards)
    [~, visVec_L] = max(cumsum(~isnan(visVec), 2), [], 2);
    refIdx = find(visVec_L==(max(visVec_L)),1);
    refShiftVal = fileRes(refIdx,shiftVal_col);
    curr_vv = visVec(refIdx,:); curr_vv = curr_vv(~isnan(curr_vv));
    
    rampV = find(abs(diff(curr_vv))~=0);
    
    % draw "ramps" & "blanks"
    for i = 1:numel(rampV)-1
        if abs(diff([diff([rampV(i),rampV(i+1)]),const.numRamp])) <= 1
            fill([rampV(i),rampV(i),rampV(i+1),rampV(i+1),rampV(i)],[0,1,1,0,0],colFade,'FaceAlpha',0.2,'EdgeAlpha',0);
        elseif diff([rampV(i),rampV(i+1)]) == const.numBlank
            fill([rampV(i),rampV(i),rampV(i+1),rampV(i+1),rampV(i)],[0,1,1,0,0],colBl,'FaceAlpha',0.2,'EdgeAlpha',0);
        end
    end
    
    % draw saccades
    for i = 1:size(fileRes,1)
  
        nbf_sacON  = floor(fileRes(i,sac_relON_col));
        nbf_sacOFF = floor(fileRes(i,sac_relOFF_col));
        
        if fileRes(i,cond_col) ~= 4
            if fileRes(i,blankType_col) == 2 % highlight most interesting category
                plot([nbf_sacON,nbf_sacOFF]-refShiftVal + fileRes(i,shiftVal_col),[1,1]*rand,'Color',colEcc{fileRes(i,ecc_col)},'LineWidth',1.5);
            else
                plot([nbf_sacON,nbf_sacOFF]-refShiftVal + fileRes(i,shiftVal_col),[1,1]*rand,'Color',[0,0,0,0.3]);
            end
        end
    end
    
    xlabel('Time in trial (s)');
    ylabel('Target visibility');
    
    set(gca,'xLim',[0,size(curr_vv,2)],'xTick',[0:24:size(visVec,2)],'xTickLabel',[0:0.2:size(visVec,2)/120]);

    
    
    %-----------------------------
    % [B] saccade duration per ecc
    %-----------------------------
    
    subplot(3,2,3); hold on;
    title(sprintf('SD: %i,  Dur: %i,  merge: %i',eyeana.velSD, eyeana.minDur, eyeana.mergeInt));
    
    for i = 1:3
        myH{i} = histogram(fileRes(fileRes(:,ecc_col)==i & fileRes(:,cond_col)~=4,sacDur_col),[0:2.5:80],'FaceColor',colEcc{i});
        maxVals(i) = max(myH{i}.Values);
        text( median(fileRes(fileRes(:,ecc_col)==i & fileRes(:,cond_col)~=4,sacDur_col)),max(myH{i}.Values)+2,sprintf('%2.0f',median(fileRes(fileRes(:,ecc_col)==i & fileRes(:,cond_col)~=4,sacDur_col))),'HorizontalAlignment','center','Color',colEcc{i});
    end
    
    ylim([0,max(maxVals)*1.25]);
    xlabel('Saccade duration (ms)');
    ylabel('Freq (abs)');
    
    
    %-----------------------
    % blank category per ecc
    %-----------------------
    subplot(3,2,4); hold on;
    title(sprintf('Ecc: %i, %i, %i dva',const.ecc_deg));
 
    for i = 1:3
        myH2{i} = histogram(fileRes(fileRes(:,ecc_col)==i & fileRes(:,cond_col)~=4,blankType_col)+(i-1)*4,'FaceColor',colEcc{i});
        maxVals2(i) = max(myH2{i}.Values);
    end
    
    xlim([-1,14]);
    set(gca,'xLim',[-1,14],'xTick',[1:12],'xTickLabel',{'11','10*','00','01'},'xTickLabelRotation',55);
    xlabel('Blank Type');
    ylabel('Freq (abs)');
    
    
    %--------------------
    % sacON re next blank
    %--------------------
    subplot(3,2,6); hold on;
    
 
    for i = 1:3
        myH3{i} = histogram(fileRes(fileRes(:,ecc_col)==i & fileRes(:,blank_tt_col)==2 & fileRes(:,cond_col)~=4,sacON_reBl_col),[-20:1:0],'FaceColor',colEcc{i});
        maxVals3(i) = max(myH3{i}.Values);
        text( mean(fileRes(fileRes(:,ecc_col)==i & fileRes(:,blank_tt_col)==2 & fileRes(:,cond_col)~=4,sacON_reBl_col)),max(myH3{i}.Values)+4,sprintf('%2.0f',median(fileRes(fileRes(:,ecc_col)==i & fileRes(:,blank_tt_col)==2 & fileRes(:,cond_col)~=4,sacON_reBl_col))),'HorizontalAlignment','center','Color',colEcc{i});
    end
    
    ylim([0,max(maxVals3)*1.25]);
    xlabel('sacON re next blank (nbf)');
    ylabel('Freq (abs)');
    
    saveas(fig0,sprintf('../../../Figures/fig0_%s.png',sub.ini));
    
    
end
    
    
%% Extraction 1 : blank duration & displ
% col#01 = split1: gap size
% col#02 = split2: displacement 
% col#03 = split3: NaN
% col#04 = split4: NaN
% col#05 = split5: NaN
% col#06 = mean response 
% col#07 = saccade latency
% col#08 = saccade amplitude
% col#09 = saccade landing error
% col#10 = number of trials

exType = '1_gap_displ';
dataFile = fileRes;

all.varTabR=[];     all.mean_resp=[];
all.mean_sacLat=[]; all.mean_sacAmp=[]; all.mean_sacErr=[];
all.num=[];         


for gap = 1:6 % [ fm07-1 ben: 1:6?]
    index.gap = dataFile(:,blankType_col) == gap;
    for displ = unique(fileRes(:,displ2_col))'
        index.displ = dataFile(:,displ2_col) == displ;
        
        allData      =   dataFile(index.gap & index.displ,:);
        val_resp     =   allData(:,resp2_col);
        val_sacLat   =   allData(:,sacLat_col);
        val_sacAmp   =   allData(:,sacAmp_col);
        val_sacErr   =   allData(:,sacErr_col);
        
        if isempty(val_resp)
            mean_resp     =   NaN;
            mean_sacLat   =   NaN;
            mean_sacAmp   =   NaN;
            mean_sacErr   =   NaN;
            num           =   0;
            varTab        =   [gap,displ];
            varTabR       =   varTab;
        else
            mean_resp     =   mean(val_resp);
            mean_sacLat   =   median(val_sacLat);
            mean_sacAmp   =   median(val_sacAmp);
            mean_sacErr   =   median(val_sacErr);
            num           =   numel(val_resp);
            varTab        =   [allData(:,blankType_col),allData(:,displ2_col)];
            varTabR       =   mean(varTab,1);
        end
        
        all.varTabR       =   [all.varTabR;       varTabR      ];
        all.mean_resp     =   [all.mean_resp;     mean_resp    ];
        all.mean_sacLat   =   [all.mean_sacLat;   mean_sacLat  ];
        all.mean_sacAmp   =   [all.mean_sacAmp;   mean_sacAmp  ];
        all.mean_sacErr   =   [all.mean_sacErr;   mean_sacErr  ];
        all.num           =   [all.num;           num          ];
    end
end


tabRes1 = [ nan(size(all.varTabR,1),5), ...
            all.mean_resp,...
            all.mean_sacLat, all.mean_sacAmp, all.mean_sacErr,...
            all.num];
tabRes1(1:size(all.varTabR,1),1:size(all.varTabR,2)) = all.varTabR;

csvwrite(sprintf('../Extract/%s_Res%s.csv',sub.ini,exType),tabRes1);




%% Plot perceptual reports
% fraction of inward / outward reports per blanking condition

if plot_1
    
    %-%-%-%-%-%-%
    
    % colors
    colAll      = [0.3,0.3,0.3];
    colPreGap   = [170, 120,  57]./255;
    colGap      = [ 41,  82, 109]./255;
    colBlank    = [100, 77,  111]./255;
    colNoBlank  = [111, 23,  7]./255;
    % [ fm07-1 ]
    
    scalfF      = 5; % scaling factor for marker size to visualize trial number (arbitrary number)
    my_fontS    = 9; % font size
    
    condTxt = {'BG','TAR','STRIPES',[],'PATCH','LANDMARK'};
    
    % read out x/y data
    xVal        = unique(fileRes(:,displ_col))';
    preGapData	= tabRes1(tabRes1(:,1)==1,6)';      preGapNum = tabRes1(tabRes1(:,1)==1,end)';
    gapData  	= tabRes1(tabRes1(:,1)==2,6)';      gapNum = tabRes1(tabRes1(:,1)==2,end)';
    BlankData  	= tabRes1(tabRes1(:,1)==6,6)';      BlankNum = tabRes1(tabRes1(:,1)==6,end)';
    noBlankData = tabRes1(tabRes1(:,1)==5,6)';      noBlankNum = tabRes1(tabRes1(:,1)==5,end)';
    % [ fm07-1 ]
    
    
    %-%-%-%-%-%-%


    % open figure
    fig1 = figure; hold on;
    title(sprintf('%s - %s (%i blocks)',const.sjct_name,condTxt{const.condition},max(fileRes(:,block_col))));
    
    
    % plot 'background lines'
    plot([0,0],[0,1],'-','Color',colAll); % midline
    plot([-100,100]/const.ppd,[0.5,0.5],'--','Color',colAll); % PSE (50%)
       
    % plot data (lines) 
    myLine_preGap   = plot(xVal,preGapData,'-','Color',colPreGap,'LineWidth',2);
    myLine_gap      = plot(xVal,gapData,   '-','Color',colGap,   'LineWidth',2);
    myLine_Blank    = plot(xVal,BlankData, '-','Color',colBlank,   'LineWidth',2);
    myLine_noBlank  = plot(xVal,noBlankData, '-','Color',colNoBlank,   'LineWidth',2);
    % [ fm07-1 ]
    
    % add markers (size represents trial number)
    scatter(xVal,preGapData,preGapNum*scalfF,'MarkerEdgeColor',colPreGap,'MarkerFaceColor',[1,1,1]);
    scatter(xVal,gapData,   gapNum*scalfF,   'MarkerEdgeColor',colGap,   'MarkerFaceColor',[1,1,1]);
    scatter(xVal,BlankData, BlankNum*scalfF, 'MarkerEdgeColor',colBlank,   'MarkerFaceColor',[1,1,1]);
    scatter(xVal,noBlankData, noBlankNum*scalfF, 'MarkerEdgeColor',colNoBlank,   'MarkerFaceColor',[1,1,1]);
    % [ fm07-1 ]
    
    
    % just for esthetics 
    scatter(xVal,preGapData,preGapNum*scalfF,'MarkerEdgeColor',colPreGap,'MarkerFaceColor',colPreGap,'MarkerFaceAlpha',0.2);
    scatter(xVal,gapData,   gapNum*scalfF,   'MarkerEdgeColor',colGap,   'MarkerFaceColor',colGap,   'MarkerFaceAlpha',0.2);
    scatter(xVal,BlankData,   BlankNum*scalfF, 'MarkerEdgeColor',colBlank, 'MarkerFaceColor',colBlank, 'MarkerFaceAlpha',0.2);
    scatter(xVal,noBlankData,   noBlankNum*scalfF,   'MarkerEdgeColor',colNoBlank,   'MarkerFaceColor',colNoBlank,   'MarkerFaceAlpha',0.2);
    % [ fm07-1 ]
    
    
    %-%-%-%-%-%-%

        
    % figure legend
    legend([myLine_preGap,myLine_gap,myLine_Blank,myLine_noBlank],...
            sprintf('preGap (%1.1f tpdp)',mean(preGapNum)),...
            sprintf('gap (%1.1f tpdp)',   mean(gapNum)),...
            sprintf('Blank (%1.1f tpdp)',   mean(BlankNum)),...
            sprintf('noBlank (%1.1f tpdp)',   mean(noBlankNum)),...
            'Location','southeast','FontSize',my_fontS);
    
    % labels    
    xlabel('INWARD   -   displacement (dva)   -   OUTWARD');
    ylabel('Report OUTWARD');
    set(gca,'YLim',[0,1],'yTick',[0:0.25:1],'xTick',xVal);
    
    % save
    saveas(fig1,sprintf('../../../Figures/fig1fm07_%s.png',sub.ini));
    
end

end