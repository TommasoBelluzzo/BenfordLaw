% [INPUT]
% data = A numeric array representing the sample on which the Benford's Law analysis must be performed.
% dran = A string representing the range of values to consider (optional, default='ALL').
%        Its value can be one of the following:
%         - ALL (all values)
%         - NEG (only negative values)
%         - POS (only positive values)
% ddec = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
%        No rounding is performed, the exceeding decimals are truncated as if they were not present.
% a    = A float [0.01,0.10] representing the statistical significance threshold for the tests (optional, default=0.05).
% ccf  = A boolean indicating whether to apply a continuity correction factor on Z-scores (optional, default=true).

function benford_analyse(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('dran','ALL',@(x)any(validatestring(x,{'ALL','NEG','POS'})));
        p.addOptional('ddec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
        p.addOptional('ccf',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    dran = res.dran;
    ddec = res.ddec;
    a = res.a;
    ccf = res.ccf;

    data_orig = data;
    data = benford_data(data,dran,ddec);

    benford_analyse_internal(data_orig,data,dran,ddec,a,ccf);

end

function benford_analyse_internal(data_orig,data,dran,ddec,a,ccf)

    [tab_1st,test_1st] = benford_extract(data,dran,ddec,'1ST',a,ccf,false);
%     gof_1st = benford_gof_all(res_1st,a);
%     
    [tab_2nd,test_2nd] = benford_extract(data,dran,ddec,'2ND',a,ccf,false);
%     gof_2nd = benford_gof_all(res_2nd,a);
%     
%     res_3rd = benford_extract(data,dran,ddec,'3RD',a,ccf,false);
%     gof_3rd = benford_gof_all(res_3rd,a);
% 
%     res_f2d = benford_extract(data,dran,ddec,'F2D',a,ccf,false);
%     gof_f2d = benford_gof_all(res_f2d,a);
% 
%     res_f3d = benford_extract(data,dran,ddec,'F3D',a,ccf,false);
%     gof_f3d = benford_gof_all(res_f3d,a);
%     
%     res_l2d = benford_extract(data,dran,ddec,'L2D',a,ccf,false);
%     gof_l2d = benford_gof_all(res_l2d,a);

	[d10n,d0n,d0p,d10p,z] = benford_duplication(data_orig,ddec);
    plot_duplication(d10n,d0n,d0p,d10p,z,ddec);

    [mant,mant_test,mant_desc] = benford_mantissae(data,dran,ddec,a,false);
    plot_mantissae(mant,mant_test,mant_desc,a);

	[df,df_coll] = benford_df(data_orig,ddec,a);
    plot_df(df,df_coll,a);

end

function plot_df(df,coll,a)

    if (df.Significance)
        df_sig = sprintf('Significant, p-Value: %.4f',df.pValue);
    else
        df_sig = sprintf('Non-Significant, p-Value: %.4f',df.pValue);
    end

    if (df.Value > 0)
        df_res = ['DF Test: ' sprintf('%.2f%%',(df.Value * 100)) ' Overstatement (' df_sig ')'];
    elseif (df.Value < 0)
        df_res = ['DF Test: ' sprintf('%.2f%%',(df.Value * -100)) ' Understatement (' df_sig ')'];
    else
        df_res = 'DF Test: No Distortion';
    end

    fig = figure();
    set(fig,'Name','Distortion Factor Model','Units','normalized','Position',[100 100 0.75 0.75]);

    sub_1 = subplot(13,9,[10 104.5]);
    histogram(sub_1,coll,10:10:100,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659],'Normalization','probability');
    set(sub_1,'XLim',[10 100]);
    title(sub_1,'Collapsed Values');
    
    sub_2 = subplot(13,9,[16 54],'Box','on');
    line([0 1],repmat(39.0865,1,2),'Color','r');
    hold on;
        line([0 1],repmat(df.EM,1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
        l1 = area(1,NaN,'FaceColor','r');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659]);
    hold off;
    set(sub_2,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_2,'YAxisLocation','right','YLim',[0 100],'YTick',[0 25 50 75 100],'YTickLabel',sprintfc('%.0f',[0 25 50 75 100].'));
    title(sub_2,sprintf('EM: %.4f | EM Max: 39.0865',df.EM));

    sub_3 = subplot(13,9,[70 108],'Box','on');
    line([0 1],repmat(df.EM,1,2),'Color','r');
    hold on;
        line([0 1],repmat(df.AM,1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
    hold off;
    set(sub_3,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_3,'YAxisLocation','right','YLim',[0 100],'YTick',[0 25 50 75 100],'YTickLabel',sprintfc('%.0f',[0 25 50 75 100].'));
    title(sub_3,sprintf('AM: %.4f | EM: %.4f',df.AM,df.EM));

    l = legend([l1 l2],'Theorical Values','Empirical Values','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    suptitle(sprintf('Distortion Factor Model (a=%.2f)\n%s',a,df_res));
    movegui(fig,'center');

end

function plot_duplication(d10n,d0n,d0p,d10p,z,ddec)

    d10n_h = height(d10n);
    d0n_h = height(d0n);
    d0p_h = height(d0p);
    d10p_h = height(d10p);
    dup = d10n_h + d0n_h + d0p_h + d10p_h;
    fs = ['%.' num2str(ddec) 'f'];

    fig = figure();
    set(fig,'Name','Number Duplication Analysis','Units','normalized','Position',[100 100 0.75 0.75]);

    if (d10p_h == 0)
        sub_1 = subplot(13,9,[19 58],'Box','on');
        line([0 1],[0 1],'Color','k','LineWidth',1.5);
        hold on;
            line([0 1],[1 0],'Color','k','LineWidth',1.5);
        hold off;
        set(sub_1,'Color',[0.941 0.941 0.941]);
        set(sub_1,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
        set(sub_1,'YLim',[0 1],'YTick',[],'YTickLabel',[]);
        title(sub_1,sprintf('Top 10 (X > 10) | No Duplicates'));
    else
        d10p_m = min([10 d10p_h]);
        d10p_d = 10 - d10p_m;
        c = [d10p.Count(1:d10p_m); NaN(d10p_d,1)];
        v = [arrayfun(@(x)num2str(x,fs),d10p.Value(1:d10p_m),'UniformOutput',false); repmat({'N/A'},d10p_d,1)];
        
        sub_1 = subplot(13,9,[19 58]);
        barh(flipud(c),'FaceColor',[0.239 0.149 0.659]);
        set(sub_1,'YLim',[0 11],'YTickLabel',flipud(v));
        title(sub_1,sprintf('Top 10 (X > 10) | Duplicates: %d',d10p_h));
    end

    if (d0p_h == 0)
        sub_2 = subplot(13,9,[24 63],'Box','on');
        line([0 1],[0 1],'Color','k','LineWidth',1.5);
        hold on;
            line([0 1],[1 0],'Color','k','LineWidth',1.5);
        hold off;
        set(sub_2,'Color',[0.941 0.941 0.941]);
        set(sub_2,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
        set(sub_2,'YLim',[0 1],'YTick',[],'YTickLabel',[]);
        title(sub_2,sprintf('Top 10 (10 > X > 0) | No Duplicates'));
    else
        d0p_m = min([10 d0p_h]);
        d0p_d = 10 - d0p_m;
        c = [d0p.Count(1:d0p_m); NaN(d0p_d,1)];
        v = [arrayfun(@(x)num2str(x,fs),d0p.Value(1:d0p_m),'UniformOutput',false); repmat({'N/A'},d0p_d,1)];
        
        sub_2 = subplot(13,9,[24 63]);
        barh(flipud(c),'FaceColor',[0.239 0.149 0.659]);
        set(sub_2,'YLim',[0 11],'YTickLabel',flipud(v));
        title(sub_2,sprintf('Top 10 (10 > X > 0) | Duplicates: %d',d0p_h));
    end
    
    if (d10n_h == 0)
        sub_3 = subplot(13,9,[73 112],'Box','on');
        line([0 1],[0 1],'Color','k','LineWidth',1.5);
        hold on;
            line([0 1],[1 0],'Color','k','LineWidth',1.5);
        hold off;
        set(sub_3,'Color',[0.941 0.941 0.941]);
        set(sub_3,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
        set(sub_3,'YLim',[0 1],'YTick',[],'YTickLabel',[]);
        title(sub_3,sprintf('Top 10 (X < -10) | No Duplicates'));
    else
        d10n_m = min([10 d10n_h]);
        d10n_d = 10 - d10p_m;
        c = [d10n.Count(1:d10n_m); NaN(d10n_d,1)];
        v = [arrayfun(@(x)num2str(x,fs),d10n.Value(1:d10n_m),'UniformOutput',false); repmat({'N/A'},d10n_d,1)];
        
        sub_3 = subplot(13,9,[73 112]);
        barh(flipud(c),'FaceColor',[0.239 0.149 0.659]);
        set(sub_3,'YLim',[0 11],'YTickLabel',flipud(v));
        title(sub_3,sprintf('Top 10 (X < -10) | Duplicates: %d',d10n_h));
    end

    if (d0n_h == 0)
        sub_4 = subplot(13,9,[78 117],'Box','on');
        line([0 1],[0 1],'Color','k','LineWidth',1.5);
        hold on;
            line([0 1],[1 0],'Color','k','LineWidth',1.5);
        hold off;
        set(sub_4,'Color',[0.941 0.941 0.941]);
        set(sub_4,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
        set(sub_4,'YLim',[0 1],'YTick',[],'YTickLabel',[]);
        title(sub_4,sprintf('Top 10 (0 > X > -10) | No Duplicates'));
    else
        d0n_m = min([10 d0n_h]);
        d0n_d = 10 - d0n_m;
        c = [d0n.Count(1:d0n_m); NaN(d0n_d,1)];
        v = [arrayfun(@(x)num2str(x,fs),d0n.Value(1:d0n_m),'UniformOutput',false); repmat({'N/A'},d0n_d,1)];
        
        sub_4 = subplot(13,9,[78 117]);
        barh(flipud(c),'FaceColor',[0.239 0.149 0.659]);
        set(sub_4,'YLim',[0 11],'YTickLabel',flipud(v));
        title(sub_4,sprintf('Top 10 (0 > X > -10) | Duplicates: %d',d0n_h));
    end

    st = suptitle(sprintf('Number Duplication Analysis\nTotal Duplicates: %d\nTotal Zeros: %d',dup,z));
    st_pos = get(st,'Position');
    set(st,'Position',[st_pos(1) -0.08 st_pos(3)]);
    
    movegui(fig,'center');

end

function plot_mantissae(mant,test,desc,a)

    mant_len = height(mant);
    the = [0; cumsum(ones(mant_len,1) ./ mant_len)];
    emp = [0; sort(mant.Mantissae)];

    [hist_cnt,hist_edg] = histcounts(mant.Mantissae,25,'BinLimits',[0,1]);
    hist_cent = (hist_edg(1:end-1) + hist_edg(2:end)) ./ 2;
    hist_avg = sum(hist_cnt) / 25;
    
    max_skew = (floor(max(abs(desc{3,:}))) + 1);
    max_kurt = (floor(max(abs(desc{4,:}))) + 1);

    if (test.H0)
        test_res = 'H0';
    else
        test_res = 'H1';
    end

    fig = figure();
    set(fig,'Name','Mantissae Analysis','Units','normalized','Position',[100 100 0.75 0.75]);

    sub_1 = subplot(13,9,[10 40.5]);
    plot(0,0,'Color','r','LineStyle','none','Marker','x','MarkerSize',15);
    hold on;
        plot(mant.X,mant.Y,'Color',[0.239 0.149 0.659],'LineStyle','none','Marker','.','MarkerSize',6);
        plot(mean(mant.X),mean(mant.Y),'Color',[0.239 0.149 0.659],'LineStyle','none','Marker','.','MarkerSize',15);
    hold off;
    set(sub_1,'XLim',[-1.1 1.1],'XTick',[-1 -0.5 0 0.5 1]);
    set(sub_1,'YLim',[-1.1 1.1],'YTick',[-1 -0.5 0 0.5 1]);
    axis equal;
    grid on;
    title(sub_1,'Arc');
    
    sub_2 = subplot(13,9,[15 45]);
    plot(the,'-r');
    hold on;
        plot(emp,'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
        l1 = area(1,NaN,'FaceColor','r');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659]);
    hold off;
    set(sub_2,'XLim',[0 mant_len]);
    set(sub_2,'YLim',[0 1]);
    title(sub_2,'Distribution');

    sub_3 = subplot(13,9,[55 103.5]);
    bar(hist_cent,hist_cnt,'hist');
    hold on;
        line([0 1],[hist_avg hist_avg],'Color','r','LineWidth',1.5);
    hold off;
    set(sub_3,'XLim',[0 1]);
    title(sub_3,'Histogram');

    sub_4 = subplot(13,9,[69 79]);
    line([0 1],repmat(desc{1,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{1,2},1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
    hold off;
    set(sub_4,'Box','on');
    set(sub_4,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_4,'YLim',[0 1],'YTick',[0 1],'YTickLabel',sprintfc('%.1f',[0 1].'));
    title(sub_4,sprintf('Mean: %.4f\nExpected: %.4f',desc{1,2},desc{1,1}));
    
    sub_5 = subplot(13,9,[71 81]);
    set(sub_5,'Box','on');
    line([0 1],repmat(desc{2,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{2,2},1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
    hold off;
    set(sub_5,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_5,'YAxisLocation','right','YLim',[0 0.25],'YTick',[0 0.25],'YTickLabel',sprintfc('%.2f',[0 0.25].'));
    title(sub_5,sprintf('Variance: %.4f\nExpected: %.4f',desc{2,2},desc{2,1}));

    sub_6 = subplot(13,9,[96 106]);
    set(sub_6,'Box','on');
    line([0 1],repmat(desc{3,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{3,2},1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
    hold off;
    set(sub_6,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_6,'YLim',[-max_skew max_skew],'YTick',[-max_skew 0 max_skew],'YTickLabel',sprintfc('%d',[-max_skew 0 max_skew].'));
    title(sub_6,sprintf('Skewness: %.4f\nExpected: %.4f',desc{3,2},desc{3,1}));
    
    sub_7 = subplot(13,9,[98 108]);
    set(sub_7,'Box','on');
    line([0 1],repmat(desc{4,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{4,2},1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
    hold off;
    set(sub_7,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_7,'YAxisLocation','right','YLim',[-max_kurt max_kurt],'YTick',[-max_kurt 0 max_kurt],'YTickLabel',sprintfc('%d',[-max_kurt 0 max_kurt].'));
    title(sub_7,sprintf('Kurtosis: %.4f\nExpected: %.4f',desc{4,2},desc{4,1}));

    l = legend([l1 l2],'Theorical Values','Empirical Values','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);
    
    suptitle(sprintf('Mantissae Analysis (a=%.2f)\nArc Test: %s (Statistic: %.4f | p-Value: %.4f)',a,test_res,test.Statistic,test.pValue));
    movegui(fig,'center');

end