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

    bar = waitbar(0,'The dataset is being analysed...');

    try
        [tab_1st,mad_1st,ssd_1st] = benford_extract(data,dran,ddec,'1ST',a,ccf,false);
        gof_1st = benford_gof_all(tab_1st,a);

        [tab_2nd,mad_2nd,ssd_2nd] = benford_extract(data,dran,ddec,'2ND',a,ccf,false);
        gof_2nd = benford_gof_all(tab_2nd,a);

        [tab_3rd,mad_3rd,ssd_3rd] = benford_extract(data,dran,ddec,'3RD',a,ccf,false);
        gof_3rd = benford_gof_all(tab_3rd,a);

        [tab_f2d,mad_f2d,ssd_f2d] = benford_extract(data,dran,ddec,'F2D',a,ccf,false);
        gof_f2d = benford_gof_all(tab_f2d,a);

        [tab_l2d,mad_l2d,ssd_l2d] = benford_extract(data,dran,ddec,'L2D',a,ccf,false);
        gof_l2d = benford_gof_all(tab_l2d,a);

        tab_so = benford_second_order(data,dran,ddec,false);
        tab_sum = benford_summation(data_orig,ddec);

        [d10n,d0n,d0p,d10p,z] = benford_duplication(data_orig,ddec);
        [mant,mant_test,mant_desc] = benford_mantissae(data,dran,ddec,a,false);
        [df,df_coll] = benford_df(data_orig,ddec,a);
        [zipf_tab,zipf_fit] = benford_zipf(data,dran,ddec,false);
        
        waitbar(100,bar,'Analysis completed!');
        pause(1);
        delete(bar);
    catch e
        delete(bar);
        rethrow(e);
    end
    
    plot_overview(data_orig,ddec);

    plot_frequencies('First Digits',tab_1st,a);
    plot_conformity('First Digits',tab_1st,mad_1st,ssd_1st,gof_1st,a);

    if (tab_2nd.EmpP(1) ~= 1)
        plot_frequencies('Second Digits',tab_2nd,a);
        plot_conformity('Second Digits',tab_2nd,mad_2nd,ssd_2nd,gof_2nd,a);
    end

    if (tab_3rd.EmpP(1) ~= 1)
        plot_frequencies('Third Digits',tab_3rd,a);
        plot_conformity('Third Digits',tab_3rd,mad_3rd,ssd_3rd,gof_3rd,a);
    end

    plot_frequencies('First-Two Digits',tab_f2d,a);
    plot_conformity('First-Two Digits',tab_f2d,mad_f2d,ssd_f2d,gof_f2d,a);

    if (tab_l2d.EmpP(1) ~= 1)
        plot_frequencies('Last-Two Digits',tab_l2d,a);
        plot_conformity('Last-Two Digits',tab_l2d,mad_l2d,ssd_l2d,gof_l2d,a);
    end
    
    plot_second_order(tab_so);
    plot_summation(tab_sum,ddec);

    plot_duplication(d10n,d0n,d0p,d10p,z,ddec);
    plot_mantissae(mant,mant_test,mant_desc,a);
    plot_df(df,df_coll,a);
    plot_zipf(zipf_tab,zipf_fit);

end

function t = figure_title(str)

    fig = gcf();
    fig_fts = get(fig,'DefaultAxesFontSize') + 4;
    fig_uni = get(fig,'Units');
    
    if (~strcmp(fig_uni,'pixels'))
        set(fig,'Units','pixels');
        fig_pos = get(fig,'Position');
        set(fig,'Units',fig_uni);
    else
        fig_pos = get(fig,'Position');
    end

    ff = ((fig_fts - 4) * 6.35) / fig_pos(4);

    tit = NaN;
    y_max = 0;
    y_min = 1;

    h = findobj(fig,'Type','axes');
    h_len = length(h);
    h_pos = zeros(h_len,4);

    for i = 1:h_len
        h_cur = h(i);
        
        fig_pos = get(h_cur,'Position');
        h_pos(i,:) = fig_pos;

        if (~strcmp(get(h_cur,'Tag'),'suptitle'))
            fig_y = fig_pos(2);
            fig_hei = fig_pos(4);
            
            if (fig_y < y_min)
                y_min = fig_y - (ff / 15);
            end

            if ((fig_hei + fig_y) > y_max)
                y_max = fig_hei + fig_y + (ff / 10);
            end
        else
            tit = h_cur;
        end
    end

    if (y_max > 0.92)
        scl = (0.92 - y_min) / (y_max - y_min);

        for i = 1:h_len
            fig_pos = h_pos(i,:);
            fig_pos(2) = ((fig_pos(2) - y_min) * scl) + y_min;
            fig_pos(4) = (fig_pos(4) * scl) - ((1 - scl) * (ff / 15));

            set(h(i),'Position',fig_pos);
        end
    end

    np = get(fig,'NextPlot');
    set(fig,'NextPlot','add');

    if (ishghandle(tit))
        delete(tit);
    end

    axes('Position',[0 1 1 1],'Tag','suptitle','Visible','off');
    t_int = text(0.50,-0.05,str,'HorizontalAlignment','center','FontSize',fig_fts);

    set(fig,'NextPlot',np);

    axes(gca());

    if (nargout == 1)
        t = t_int;
    end

end

function plot_conformity(name,tab,mad,ssd,gof,a)

    len = height(tab);
    seq = 1:len;
    seq_plot = 0.5:(len + 0.5);

    if (strcmp(name,'Last-Two Digits'))
        labs = sprintfc('%02d',tab.Digits);
        xrot = 90;
    else
        labs = sprintfc('%d',tab.Digits);
        
        if (len <= 10)
            xrot = 0;
        else
            xrot = 90;
        end
    end

    cv = norminv(1 - (a / 2));

    z_ok = tab.Z;
    z_ok(~tab.ZTest) = NaN;
    z_ko = tab.Z;
    z_ko(tab.ZTest) = NaN;
    
    z_fail_sum = sum(~tab.ZTest);
    z_fail_prc = (z_fail_sum / len) * 100;
    
    if (z_fail_sum == 1)
        z_fail_char = '';
    else
        z_fail_char = 's';
    end

    z_max = max(tab.Z);
    ymax = max([z_max (ceil(cv * 10) / 10)]);
    
    gof_fail = ~gof.H0;
    gof_fail_sum = sum(gof_fail);
    gof_len = height(gof);

    gof_ok = ones(gof_len,1);
    gof_ok(gof_fail) = 0;
    gof_ko = zeros(gof_len,1);
    gof_ko(gof_fail) = 1;

    mad_res = mad(~isnan(mad.Value),{'Conformity' 'Value'});
    mad_thr = arrayfun(@(x)sprintf('%.4f (%s)',mad.Threshold(x),regexprep(mad.Conformity{x},'[a-z ]+','')),1:3,'UniformOutput',false);
    
    if (isnan(mad.Value(4)))
        mad_clr = [0.239 0.149 0.659];
    else
        mad_clr = [1.000 0.840 0.000];
    end
    
    ssd_res = ssd(~isnan(ssd.Value),{'Conformity' 'Value'});
    ssd_thr = arrayfun(@(x)sprintf('%d (%s)',ssd.Threshold(x),regexprep(ssd.Conformity{x},'[a-z ]+','')),1:3,'UniformOutput',false);

    if (isnan(ssd.Value(4)))
        ssd_clr = [0.239 0.149 0.659];
    else
        ssd_clr = [1.000 0.840 0.000];
    end

    fig = figure();
    set(fig,'Name',[name ' Analysis: Conformity Analysis'],'Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 54]);
    set(bar(z_ok,'hist'),'FaceColor',[0.239 0.149 0.659]);
    hold on;
        set(bar(z_ko,'hist'),'FaceColor',[1.000 0.840 0.000]);
        line([seq_plot(1) seq_plot(end)],[cv cv],'Color','r','LineWidth',1.5);
        line([seq_plot(1) seq_plot(end)],[ymax ymax],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        line([seq_plot(end) seq_plot(end)],[0 ymax],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        l1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659],'ShowBaseLine','off');
        l3 = area(1,NaN,'FaceColor',[1.000 0.840 0.000],'ShowBaseLine','off');
    hold off;
    set(sub_1,'XLim',[seq_plot(1) seq_plot(end)],'XTick',seq,'XTickLabel',labs,'XTickLabelRotation',xrot);
    set(sub_1,'YLim',[0 ymax]);
    set(sub_1,'Box','off','TickDir','both');
    t1 = title(sub_1,sprintf('Z-scores (%d Failed Test%s, %.2f%%)',z_fail_sum,z_fail_char,z_fail_prc));
    set(t1,'Units','normalized');
    t1_pos = get(t1,'Position');
    set(t1,'Position',[0.4783 t1_pos(2) t1_pos(3)]);

    sub_2 = subplot(13,9,[64 104]);
    set(bar(gof_ok,'hist'),'FaceColor',[0.239 0.149 0.659]);
    hold on;
        set(bar(gof_ko,'hist'),'FaceColor',[1.000 0.840 0.000]);
    hold off;
    set(sub_2,'XLim',[0.5 (gof_len + 0.5)],'XTickLabel',gof.Properties.RowNames);
    set(sub_2,'YLim',[0 1],'YTick',[],'YTickLabel',[]);
    title(sub_2,sprintf('Goodness-of-Fit (Failed Tests: %d/%d)',gof_fail_sum,gof_len));
    
    sub_3 = subplot(13,9,[70 81]);
    set(bar(1,mad_res.Value,0.5),'FaceColor',mad_clr);
    hold on;
        for i = 1:3,line([0.5 1.5],[mad.Threshold(i) mad.Threshold(i)],'Color','r'),end
    hold off;
    set(sub_3,'XLim',[0.5 1.5],'XTick',[],'XTickLabel',[]);
    set(sub_3,'YLim',[0.0 (mad.Threshold(3) * 1.2)],'YTick',mad.Threshold(1:3),'YTickLabel',mad_thr);
    title(sub_3,sprintf('MAD Test: %s (%.4f)',mad_res.Conformity{:},mad_res.Value));
    
    sub_4 = subplot(13,9,[97 108]);
    set(bar(1,ssd_res.Value,0.5),'FaceColor',ssd_clr);
    hold on;
        for i = 1:3,line([0.5 1.5],[ssd.Threshold(i) ssd.Threshold(i)],'Color','r'),end
    hold off;
    set(sub_4,'XLim',[0.5 1.5],'XTick',[],'XTickLabel',[]);
    set(sub_4,'YLim',[0.0 (ssd.Threshold(3) * 1.2)],'YTick',ssd.Threshold(1:3),'YTickLabel',ssd_thr);
    title(sub_4,sprintf('SSD Test: %s (%.2f)',ssd_res.Conformity{:},ssd_res.Value));

    l = legend([l1 l2 l3],'Theorical Values','Empirical Values','Failed Tests','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    set(fig,'UserData',gof);
    dcm = datacursormode(fig); 
    set(dcm,'Enable','on'); 
    set(dcm,'UpdateFcn',@plot_conformity_tooltip);
    
    figure_title(sprintf('%s Analysis (a=%.2f)\nConformity Analysis',name,a));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function out = plot_conformity_tooltip(~,evt_obj)

    targ = get(evt_obj.Target,'Parent');
    targ_tit = get(targ,'Title');
    targ_txt = get(targ_tit,'String');

    if (~startsWith(targ_txt,'Goodness-of-Fit'))
        out = '';
        return;
    end

    fig = get(targ,'Parent');
    pos = get(evt_obj,'Position');

    gof = get(fig,'UserData');
    gof_off = ceil(pos(1));
    
    if ((gof_off < 1) || (gof_off > height(gof)))
        out = '';
        return;
    end 
    
    gof_row = gof(gof_off,:);

    out = num2str([gof_row.Statistics, gof_row.pValues].','%0.4f');
    out = [pad([gof_row.Properties.RowNames{:} ':'],numel(out(1,:))); out];

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
    set(fig,'Name','Distortion Factor Model','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 104.5]);
    histogram(sub_1,coll,10:10:100,'FaceAlpha',1,'FaceColor',[0.239 0.149 0.659],'Normalization','probability');
    set(sub_1,'XLim',[10 100]);
    title(sub_1,'Collapsed Values');

    sub_2 = subplot(13,9,[16 54]);
    line([0.5 1.5],repmat(39.0865,1,2),'Color','r');
    hold on;
        set(bar(1,df.EM,0.5),'FaceColor',[0.239 0.149 0.659]);
        l1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659],'ShowBaseLine','off');
    hold off;
    set(sub_2,'Box','on');
    set(sub_2,'XLim',[0.5 1.5],'XTick',[],'XTickLabel',[]);
    set(sub_2,'YLim',[0 50],'YTick',0:25:50,'YTickLabel',sprintfc('%.0f',(0:25:50).'));
    title(sub_2,sprintf('EM: %.4f | EM Max: 39.0865',df.EM));

    sub_3 = subplot(13,9,[70 108]);
    line([0.5 1.5],repmat(df.EM,1,2),'Color','r');
    hold on;
        set(bar(1,df.AM,0.5),'FaceColor',[0.239 0.149 0.659]);
    hold off;
    set(sub_3,'Box','on');
    set(sub_3,'XLim',[0.5 1.5],'XTick',[],'XTickLabel',[]);
    set(sub_3,'YLim',[0 50],'YTick',0:25:50,'YTickLabel',sprintfc('%.0f',(0:25:50).'));
    title(sub_3,sprintf('AM: %.4f | EM: %.4f',df.AM,df.EM));

    l = legend([l1 l2],'Theorical Values','Empirical Values','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    figure_title(sprintf('Distortion Factor Model (a=%.2f)\n%s',a,df_res));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_duplication(d10n,d0n,d0p,d10p,z,ddec)

    d10n_h = height(d10n);
    d0n_h = height(d0n);
    d0p_h = height(d0p);
    d10p_h = height(d10p);
    dup = d10n_h + d0n_h + d0p_h + d10p_h;
    fs = ['%.' num2str(ddec) 'f'];

    fig = figure();
    set(fig,'Name','Number Duplication Analysis','Units','normalized','Position',[100 100 0.85 0.85]);

    if (d10p_h == 0)
        sub_1 = subplot(13,9,[19 58]);
        line([0 1],[0 1],'Color','k','LineWidth',1.5);
        hold on;
            line([0 1],[1 0],'Color','k','LineWidth',1.5);
        hold off;
        set(sub_1,'Box','on','Color',[0.941 0.941 0.941]);
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
        sub_2 = subplot(13,9,[24 63]);
        line([0 1],[0 1],'Color','k','LineWidth',1.5);
        hold on;
            line([0 1],[1 0],'Color','k','LineWidth',1.5);
        hold off;
        set(sub_2,'Box','on','Color',[0.941 0.941 0.941]);
        set(sub_2,'XDir','reverse','XLim',[0 1],'XTick',[],'XTickLabel',[]);
        set(sub_2,'YAxisLocation','right','YLim',[0 1],'YTick',[],'YTickLabel',[]);
        title(sub_2,sprintf('Top 10 (10 > X > 0) | No Duplicates'));
    else
        d0p_m = min([10 d0p_h]);
        d0p_d = 10 - d0p_m;
        c = [d0p.Count(1:d0p_m); NaN(d0p_d,1)];
        v = [arrayfun(@(x)num2str(x,fs),d0p.Value(1:d0p_m),'UniformOutput',false); repmat({'N/A'},d0p_d,1)];
        
        sub_2 = subplot(13,9,[24 63]);
        barh(flipud(c),'FaceColor',[0.239 0.149 0.659]);
        set(sub_2,'XDir','reverse');
        set(sub_2,'YAxisLocation','right','YLim',[0 11],'YTickLabel',flipud(v));
        title(sub_2,sprintf('Top 10 (10 > X > 0) | Duplicates: %d',d0p_h));
    end
    
    if (d10n_h == 0)
        sub_3 = subplot(13,9,[73 112]);
        line([0 1],[0 1],'Color','k','LineWidth',1.5);
        hold on;
            line([0 1],[1 0],'Color','k','LineWidth',1.5);
        hold off;
        set(sub_3,'Box','on','Color',[0.941 0.941 0.941]);
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
        sub_4 = subplot(13,9,[78 117]);
        line([0 1],[0 1],'Color','k','LineWidth',1.5);
        hold on;
            line([0 1],[1 0],'Color','k','LineWidth',1.5);
        hold off;
        set(sub_4,'Box','on','Color',[0.941 0.941 0.941]);
        set(sub_4,'XDir','reverse','XLim',[0 1],'XTick',[],'XTickLabel',[]);
        set(sub_4,'YLim',[0 1],'YTick',[],'YTickLabel',[]);
        title(sub_4,sprintf('Top 10 (0 > X > -10) | No Duplicates'));
    else
        d0n_m = min([10 d0n_h]);
        d0n_d = 10 - d0n_m;
        c = [d0n.Count(1:d0n_m); NaN(d0n_d,1)];
        v = [arrayfun(@(x)num2str(x,fs),d0n.Value(1:d0n_m),'UniformOutput',false); repmat({'N/A'},d0n_d,1)];
        
        sub_4 = subplot(13,9,[78 117]);
        barh(flipud(c),'FaceColor',[0.239 0.149 0.659]);
        set(sub_4,'XDir','reverse');
        set(sub_4,'YAxisLocation','right','YLim',[0 11],'YTickLabel',flipud(v));
        title(sub_4,sprintf('Top 10 (0 > X > -10) | Duplicates: %d',d0n_h));
    end

    st = figure_title(sprintf('Number Duplication Analysis\nTotal Duplicates: %d\nTotal Zeros: %d',dup,z));
    st_pos = get(st,'Position');
    set(st,'Position',[0.5170 -0.080 st_pos(3)]);
    
    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_frequencies(name,tab,a)

    len = height(tab);
    seq = 1:len;
    seq_plot = 0.5:(len + 0.5);
    
    l2d = strcmp(name,'Last-Two Digits');
    
    if (l2d)
        labs = sprintfc('%02d',tab.Digits);
        xrot = 90;
    else
        labs = sprintfc('%d',tab.Digits);
        
        if (len <= 10)
            xrot = 0;
        else
            xrot = 90;
        end
    end
    
    if (strcmp(name,'First-Two Digits') || l2d)
        fd = floor(tab.Digits ./ 10);
        sd = mod(tab.Digits,10);
        diff = abs(fd - sd);
        
        tot = sum(tab.EmpC);

        adj = sum((diff == 0) | (diff == 1) | (diff == 9));
        adj_prc = (adj / tot) * 100;
        adj_txt = sprintf('Adjacent Digits: %.2f%% | ',adj_prc);
    else
        adj_txt = '';
    end

    emp_sus = ((tab.EmpP ~= 0) & (tab.EmpP < tab.ThePL)) | (tab.EmpP > tab.ThePU);
    emp_sus_prc = (sum(emp_sus) / len) * 100;

    emp_ok = tab.EmpP;
    emp_ok(emp_sus) = NaN;
    emp_ko = tab.EmpP;
    emp_ko(~emp_sus) = NaN;

    the_x = linspace(1,len,1000);
    the_pu_y = interp1(seq,tab.ThePU,the_x,'spline');
    the_p_y = interp1(seq,tab.TheP,the_x,'spline');
    the_pl_y = interp1(seq,tab.ThePL,the_x,'spline');

    fig = figure();
    set(fig,'Name',[name ' Analysis: Frequencies Analysis'],'Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 54]);
    set(bar(emp_ok,'hist'),'FaceColor',[0.239 0.149 0.659]);
    hold on;
        set(bar(emp_ko,'hist'),'FaceColor',[1.000 0.840 0.000]);
        plot(the_x,the_pu_y,'Color','r','LineStyle','--','LineWidth',1.5);
        plot(the_x,the_p_y,'Color','r','LineWidth',1.5);
        plot(the_x,the_pl_y,'Color','r','LineStyle','--','LineWidth',1.5);
        ylim = get(sub_1,'YLim');
        line([seq_plot(1) seq_plot(end)],[ylim(2) ylim(2)],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        line([seq_plot(end) seq_plot(end)],[ylim(1) ylim(2)],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        l1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659],'ShowBaseLine','off');
        l3 = area(1,NaN,'FaceColor',[1.000 0.840 0.000],'ShowBaseLine','off');
    hold off;
    set(sub_1,'XLim',[seq_plot(1) seq_plot(end)],'XTick',seq,'XTickLabel',labs,'XTickLabelRotation',xrot);
    set(sub_1,'YLim',[ylim(1) ylim(2)]);
    set(sub_1,'Box','off','TickDir','both');
    t1 = title(sub_1,sprintf('Frequencies (%sSuspicious Digits: %.2f%%)',adj_txt,emp_sus_prc));
    set(t1,'Units','normalized');
    t1_pos = get(t1,'Position');
    set(t1,'Position',[0.4783 t1_pos(2) t1_pos(3)]);
    
    sub_2 = subplot(13,9,[64 108]);
    set(bar(tab.EmpF,'hist'),'FaceColor',[0.239 0.149 0.659]);
    hold on;
        plot(tab.TheF(1:end-1),'Color','r','LineStyle','none','Marker','.','MarkerSize',15);
        for i = seq(1:end-1),line([seq_plot(i) seq_plot(i+1)],[tab.TheF(i) tab.TheF(i)],'Color','r'),end
        line([seq_plot(1) seq_plot(end)],[1 1],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        line([seq_plot(end) seq_plot(end)],[0 1],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
    hold off;
    set(sub_2,'XLim',[seq_plot(1) seq_plot(end)],'XTick',seq,'XTickLabel',labs,'XTickLabelRotation',xrot);
    set(sub_2,'YLim',[0 1]);
    set(sub_2,'Box','off','TickDir','both');
    t2 = title(sub_2,'Cumulative Frequencies');
    set(t2,'Units','normalized');
    t2_pos = get(t2,'Position');
    set(t2,'Position',[0.4783 t2_pos(2) t2_pos(3)]);

    l = legend([l1 l2 l3],'Theorical Values','Empirical Values','Suspicious Digits','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    figure_title(sprintf('%s Analysis (a=%.2f)\nFrequencies Analysis',name,a));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

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
    set(fig,'Name','Mantissae Analysis','Units','normalized','Position',[100 100 0.85 0.85]);

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
        l1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659],'ShowBaseLine','off');
    hold off;
    set(sub_2,'XLim',[0 mant_len]);
    set(sub_2,'YLim',[0 1]);
    title(sub_2,'Distribution');

    sub_3 = subplot(13,9,[55 103.5]);
    set(bar(hist_cent,hist_cnt,'hist'),'FaceColor',[0.239 0.149 0.659]);
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
    line([0 1],repmat(desc{2,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{2,2},1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
    hold off;
    set(sub_5,'Box','on');
    set(sub_5,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_5,'YAxisLocation','right','YLim',[0 0.25],'YTick',[0 0.25],'YTickLabel',sprintfc('%.2f',[0 0.25].'));
    title(sub_5,sprintf('Variance: %.4f\nExpected: %.4f',desc{2,2},desc{2,1}));

    sub_6 = subplot(13,9,[96 106]);
    line([0 1],repmat(desc{3,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{3,2},1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
    hold off;
    set(sub_6,'Box','on');
    set(sub_6,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_6,'YLim',[-max_skew max_skew],'YTick',[-max_skew 0 max_skew],'YTickLabel',sprintfc('%d',[-max_skew 0 max_skew].'));
    title(sub_6,sprintf('Skewness: %.4f\nExpected: %.4f',desc{3,2},desc{3,1}));
    
    sub_7 = subplot(13,9,[98 108]);
    line([0 1],repmat(desc{4,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{4,2},1,2),'Color',[0.239 0.149 0.659],'LineStyle',':','LineWidth',1.5);
    hold off;
    set(sub_7,'Box','on');
    set(sub_7,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_7,'YAxisLocation','right','YLim',[-max_kurt max_kurt],'YTick',[-max_kurt 0 max_kurt],'YTickLabel',sprintfc('%d',[-max_kurt 0 max_kurt].'));
    title(sub_7,sprintf('Kurtosis: %.4f\nExpected: %.4f',desc{4,2},desc{4,1}));

    l = legend([l1 l2],'Theorical Values','Empirical Values','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);
    
    figure_title(sprintf('Mantissae Analysis (a=%.2f)\nArc Test: %s (Statistic: %.4f | p-Value: %.4f)',a,test_res,test.Statistic,test.pValue));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_overview(data_orig,ddec)

    data_orig = data_orig(:);

    data_sor = sort(unique(abs(data_orig)),'descend');
    data_max = max(data_orig);
    data_min = min(data_orig);
    data_mag = floor(log10(data_sor(1)) + 1);
    data_rsf = data_sor(1) / data_sor(2);
    data_avg = mean(data_orig);
    data_std = std(data_orig);
    
    nb = round((data_max - data_min) / (2 * iqr(data_orig) * (numel(data_orig) ^ (-1/3))),0);

    fig = figure();
    set(fig,'Name','Dataset Overview','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 54]);
    set(bar(data_orig,'hist'),'FaceColor',[0.239 0.149 0.659]);
    set(sub_1,'XLim',[0.5 (numel(data_orig) + 0.5)],'XTick',[],'XTickLabel',[]);
    set(sub_1,'YLim',[data_min data_max]);
    t1 = title(sub_1,'Values');
    set(t1,'Units','normalized');
    t1_pos = get(t1,'Position');
    set(t1,'Position',[0.4783 t1_pos(2) t1_pos(3)]);

    sub_2 = subplot(13,9,[64 117]);
    h = histfit(data_orig,nb,'kernel');
    set(h(1),'FaceColor',[0.239 0.149 0.659]);
    set(h(2),'Color','r','LineWidth',0.75);
    t2 = title(sub_2,'Histogram & Kernel Density');
    set(t2,'Units','normalized');
    t2_pos = get(t2,'Position');
    set(t2,'Position',[0.4783 t2_pos(2) t2_pos(3)]);

    figure_title(sprintf(sprintf('Dataset Overview\nMaximum: %%.%df | Minimum: %%.%df | Magnitude: %%d | RSF: %%.%df | Mean: %%.%df | SD: %%.%df',ddec,ddec,ddec,ddec,ddec),data_max,data_min,data_mag,data_rsf,data_avg,data_std));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_second_order(tab)

    len = height(tab);
    seq = 1:len;
    seq_plot = 0.5:(len + 0.5);
    
    ymax = max([tab.TheP; tab.EmpP]);
    
    the_x = linspace(1,len,1000);
    the_p = interp1(seq,tab.TheP,the_x,'spline');
    
    crit = NaN(len,1);
    crit(1:10:90) = ymax;

    fig = figure();
    set(fig,'Name','Second Order Analysis','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[1 108]);
    set(bar(crit,'hist'),'FaceColor',[1.000 0.840 0.000]);
    hold on;
        set(bar(tab.EmpP,'hist'),'FaceColor',[0.239 0.149 0.659]);
        plot(the_x,the_p,'Color','r','LineWidth',1.5);
        line([seq_plot(1) seq_plot(end)],[ymax ymax],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        line([seq_plot(end) seq_plot(end)],[0 ymax],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        l1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659],'ShowBaseLine','off');
        l3 = area(1,NaN,'FaceColor',[1.000 0.840 0.000],'ShowBaseLine','off');
    hold off;
    set(sub_1,'XLim',[seq_plot(1) seq_plot(end)],'XTick',seq,'XTickLabel',sprintfc('%02d',tab.Digits),'XTickLabelRotation',90);
    set(sub_1,'YLim',[0 ymax]);
    set(sub_1,'Box','off','TickDir','both');

    l = legend([l1 l2 l3],'Theorical Values','Empirical Values','Expected Spikes','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);
    
    figure_title('Second Order Analysis');

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_summation(tab,ddec)

    len = height(tab);
    seq = 1:len;
    seq_plot = 0.5:(len + 0.5);
    
    ymax = max([tab.TheP; tab.EmpP]);

    sum_max = max(tab.Smt);
    sum_min = min(tab.Smt(tab.Smt > 0));
    avg_aes = mean(tab.AES);
    
    the_x = linspace(1,len,1000);
    the_p = interp1(seq,tab.TheP,the_x,'spline');

    fig = figure();
    set(fig,'Name','Summation Analysis','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 54]);
    set(bar(tab.EmpP,'hist'),'FaceColor',[0.239 0.149 0.659]);
    hold on;
        plot(the_x,the_p,'Color','r','LineWidth',1.5);
        line([seq_plot(1) seq_plot(end)],[ymax ymax],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        line([seq_plot(end) seq_plot(end)],[0 ymax],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        l1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659],'ShowBaseLine','off');
    hold off;
    set(sub_1,'XLim',[seq_plot(1) seq_plot(end)],'XTick',seq,'XTickLabel',sprintfc('%02d',tab.Digits),'XTickLabelRotation',90);
    set(sub_1,'YLim',[0 ymax]);
    set(sub_1,'Box','off','TickDir','both');
    t1 = title(sub_1,'Frequencies');
    set(t1,'Units','normalized');
    t1_pos = get(t1,'Position');
    set(t1,'Position',[0.4783 t1_pos(2) t1_pos(3)]);
    
    sub_2 = subplot(13,9,[64 108]);
    plot(tab.EmpF,'Color',[0.239 0.149 0.659]);
    hold on;
        plot(tab.TheF,'Color','r');
        line([1 90],[1 1],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
        line([90 90],[0 1],'AlignVertexCenters','on','Color','k','LineWidth',0.01);
    hold off;
    set(sub_2,'XLim',[1 90],'XTick',1:10:90,'XTickLabel',sprintfc('%02d',10:10:90));
    set(sub_2,'YLim',[0 1]);
    set(sub_2,'Box','off','TickDir','both');
    title(sub_2,'Cumulative Frequencies');
    
    l = legend([l1 l2],'Theorical Values','Empirical Values','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);
    
    figure_title(sprintf(sprintf('Summation Analysis\nMaximum: %%.%df | Minimum: %%.%df | Average AES: %%.%df',ddec,ddec,ddec),sum_max,sum_min,avg_aes));

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end

function plot_zipf(tab,fit)

    the_b = fit.The(2);
    
    if (sign(the_b) > 0)
        the_b_val = the_b;
        the_b_sign = '+';
    else
        the_b_val = -1 * the_b;
        the_b_sign = '-';
    end
        
    the_eqn = ['y = ' num2str(fit.The(1),'%.2f') ' ' the_b_sign ' ' num2str(the_b_val,'%.2f') 'x (R2=' num2str((fit.The(3) * 100),'%.2f%%') ')'];
    the_val = fit.The(1) + (tab.RankL .* the_b);
    
    emp_b = fit.Emp(2);

    if (sign(emp_b) > 0)
        emp_b_val = emp_b;
        emp_b_sign = '+';
    else
        emp_b_val = -1 * emp_b;
        emp_b_sign = '-';
    end

    emp_eqn = ['y = ' num2str(fit.Emp(1),'%.2f') ' ' emp_b_sign ' ' num2str(emp_b_val,'%.2f') 'x (R2=' num2str((fit.Emp(3) * 100),'%.2f%%') ')'];
    emp_val = fit.Emp(1) + (tab.RankL .* emp_b);
    
    r = linspace(0,max(tab.RankL),1000);
    the_p = interp1(tab.RankL,tab.TheP,r,'spline');
    emp_p = interp1(tab.RankL,tab.EmpP,r,'spline');
    
    xlim = max(tab.RankL);
    ylim = max([the_val; emp_val]);

    fig = figure();
    set(fig,'Name','Zipf''s Law Analysis','Units','normalized','Position',[100 100 0.85 0.85]);

    sub_1 = subplot(13,9,[10 54]);
    plot(tab.RankL,tab.SampleL,'Color',[0.239 0.149 0.659],'LineStyle','none','Marker','o');
    hold on;
        sub_l_l1 = plot(tab.RankL,the_val,'Color','r');
        sub_l_l2 = plot(tab.RankL,emp_val,'Color',[0.239 0.149 0.659]);
    hold off;
    set(sub_1,'Box','on');
    set(sub_1,'XLim',[0 xlim]);
    set(sub_1,'YLim',[0 ylim]);
    sub_1_l = legend([sub_l_l1 sub_l_l2],the_eqn,emp_eqn,'Location','best','Orientation','vertical','Units','normalized');
    set(sub_1_l,'Box','off','FontSize',10);
    t1 = title(sub_1,'Distributions');
    set(t1,'Units','normalized');
    t1_pos = get(t1,'Position');
    set(t1,'Position',[0.4783 t1_pos(2) t1_pos(3)]);

    sub_2 = subplot(13,9,[64 108]);
    plot(r,emp_p,'Color',[0.239 0.149 0.659]);
    hold on;
        plot(r,the_p,'Color','r');
        l1 = area(1,NaN,'FaceColor','r','ShowBaseLine','off');
        l2 = area(1,NaN,'FaceColor',[0.239 0.149 0.659],'ShowBaseLine','off');
    hold off;
    set(sub_2,'Box','on');
    set(sub_2,'XLim',[0 xlim]);
    set(sub_2,'YLim',[0 1],'YTick',0:0.2:1,'YTickLabel',sprintfc('%.0f%%',(0:20:100).'));
    t2 = title(sub_2,'Proportions');
    set(t2,'Units','normalized');
    t2_pos = get(t2,'Position');
    set(t2,'Position',[0.4783 t2_pos(2) t2_pos(3)]);

    l = legend([l1 l2],'Theorical Values','Empirical Values','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Box','off','Position',[((1 - l_pos(3)) / 2) 0.067 l_pos(3) l_pos(4)]);

    figure_title('Zipf''s Law Analysis');

    pause(0.01);

    jfr = get(fig,'JavaFrame');
    set(jfr,'Maximized',true);

end
