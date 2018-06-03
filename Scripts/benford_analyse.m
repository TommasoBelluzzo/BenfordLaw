% [INPUT]
% data = A numeric array representing the sample on which the Benford's Law analysis must be performed.
% ran  = A string representing the range of values to consider (optional, default='ALL').
%        Its value can be one of the following:
%         - ALL (all values)
%         - NEG (only negative values)
%         - POS (only positive values)
% dec  = An integer [0,10] representing the number of decimal places to consider (optional, default=2).
% a    = A float [0.01,0.10] representing the statistical significance threshold for the tests (optional, default=0.05).
% ccf  = A boolean indicating whether to apply a continuity correction factor on Z-scores (optional, default=true).

function benford_analyse(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('ran','ALL',@(x)any(validatestring(x,{'ALL','NEG','POS'})));
        p.addOptional('dec',2,@(x)validateattributes(x,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',10}));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
        p.addOptional('ccf',true,@(x)validateattributes(x,{'logical'},{'scalar'}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    ran = res.ran;
    dec = res.dec;
    a = res.a;
    ccf = res.ccf;

    data = benford_data(data,ran,dec);

    benford_analyse_internal(data,ran,dec,a,ccf);

end

function benford_analyse_internal(data,ran,dec,a,ccf)

    res = benford_digits(data,'L2D',a,ccf,ran,dec,false);

    [mant,mant_test,mant_desc] = benford_mantissae(data,a,ran,dec,false);
    plot_mantissae(mant,mant_test,mant_desc);

%     bd2 = benford_data(data,2);
%     gofs_fo = benford_gof_all(bd2,a,false);
%     gofs_so = benford_gof_all(bd2,a,true);

end

function plot_mantissae(mant,test,desc)

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
        line([0, 1],[hist_avg hist_avg],'Color','r','LineWidth',1.5);
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
    
    suptitle(sprintf('Mantissae Analysis\nArc Test: %s (Statistic: %.4f | p-Value: %.4f)',test_res,test.Statistic,test.pValue));
    movegui(fig,'center');

end