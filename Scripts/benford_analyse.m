% [INPUT]
% data = A numeric array representing the sample on which the Benford's Law analysis must be performed.
% a    = A float [0.01,0.10] representing the statistical significance threshold for the tests (optional, default=0.05).

function benford_analyse(varargin)

    persistent p;

    if (isempty(p))
        p = inputParser();
        p.addRequired('data',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
        p.addOptional('a',0.05,@(x)validateattributes(x,{'double','single'},{'scalar','real','finite','>=',0.01,'<=',0.10}));
    end

    p.parse(varargin{:});

    res = p.Results;
    data = res.data;
    a = res.a;
    
    benford_analyse_internal(data,a);

end

function benford_analyse_internal(data,a)


    bd2 = benford_data(data,2);

    [mant_test,mant_desc] = benford_mantissae(bd2,a);
    plot_mantissae(bd2.Mantissae,mant_test,mant_desc);
    
    gofs_fo = benford_gof_all(bd2,a,false);
    gofs_so = benford_gof_all(bd2,a,true);

end

function plot_mantissae(mant,test,desc)

    mant_len = numel(mant);
    the = cumsum(ones(mant_len,1) ./ mant_len);
    emp = sort(mant);
    
    if (test.H0)
        test_res = 'H0 Accepted';
    else
        test_res = 'H0 Rejected';
    end
    
    max_skew = (floor(max(abs(desc{3,:}))) + 1);
    max_kurt = (floor(max(abs(desc{4,:}))) + 1);

    [hist_cnt,hist_edg] = histcounts(mant,25);
    hist_cent = (hist_edg(1:end-1) + hist_edg(2:end)) ./ 2;
    hist_avg = sum(hist_cnt) / 25;

    fig = figure();
    set(fig,'Name','Mantissae Analysis','Units','normalized','Position',[100 100 0.75 0.75]);

    sub_1 = subplot(5,4,[1 10]);
    plot(the,'-r');
    hold on;
        plot(emp,'Color',[0.239 0.149 0.659],'LineStyle','--');
        l1 = area(1,NaN,'FaceColor','r');
        l2 = area(2,NaN,'FaceColor',[0.239 0.149 0.659]);
    hold off;
    set(sub_1,'XLim',[0 mant_len]);
    set(sub_1,'YLim',[0 1]);
    title(sub_1,'Distribution');
    
    sub_2 = subplot(5,4,[3 12]);
    bar(hist_cent,hist_cnt,'hist');
    hold on;
        line([0, 1],[hist_avg hist_avg],'Color','r');
    hold off;
    set(sub_2,'XLim',[0 1]);
    title(sub_2,'Histogram');
    
    sub_3 = subplot(5,4,17);
    line([0 1],repmat(desc{1,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{1,2},1,2),'Color',[0.239 0.149 0.659]);
    hold off;
    set(sub_3,'Box','on');
    set(sub_3,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_3,'YLim',[0 1],'YTick',[0 0.5 1],'YTickLabel',sprintfc('%.1f',[0 0.5 1].'));
    title(sub_3,sprintf('Mean: %.4f\nExpected: %.4f',desc{1,2},desc{1,1}));
    
    sub_4 = subplot(5,4,18);
    set(sub_4,'Box','on');
    line([0 1],repmat(desc{2,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{2,2},1,2),'Color',[0.239 0.149 0.659]);
    hold off;
    set(sub_4,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_4,'YLim',[0 0.25],'YTick',[0 0.25],'YTickLabel',sprintfc('%.2f',[0 0.25].'));
    title(sub_4,sprintf('Variance: %.4f\nExpected: %.4f',desc{2,2},desc{2,1}));

    sub_5 = subplot(5,4,19);
    set(sub_5,'Box','on');
    line([0 1],repmat(desc{3,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{3,2},1,2),'Color',[0.239 0.149 0.659]);
    hold off;
    set(sub_5,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_5,'YLim',[-max_skew max_skew],'YTick',[-max_skew 0 max_skew],'YTickLabel',sprintfc('%d',[-max_skew 0 max_skew].'));
    title(sub_5,sprintf('Skewness: %.4f\nExpected: %.4f',desc{3,2},desc{3,1}));
    
    sub_6 = subplot(5,4,20);
    set(sub_6,'Box','on');
    line([0 1],repmat(desc{4,1},1,2),'Color','r');
    hold on;
        line([0 1],repmat(desc{4,2},1,2),'Color',[0.239 0.149 0.659]);
    hold off;
    set(sub_6,'XLim',[0 1],'XTick',[],'XTickLabel',[]);
    set(sub_6,'YLim',[-max_kurt max_kurt],'YTick',[-max_kurt 0 max_kurt],'YTickLabel',sprintfc('%d',[-max_kurt 0 max_kurt].'));
    title(sub_6,sprintf('Kurtosis: %.4f\nExpected: %.4f',desc{4,2},desc{4,1}));

    l = legend([l1 l2],'Theorical Values','Empirical Values','Location','best','Orientation','horizontal','Units','normalized');
    l_pos = get(l,'Position');
    set(l,'Position',[((1 - l_pos(3)) / 2) 0.32 l_pos(3) l_pos(4)]);
    
    suptitle(sprintf('Mantissae Analysis\nMantissae Arc Test: %s | Statistic: %.4f | p-Value: %.4f',test_res,test.Statistic,test.pValue));
    movegui(fig,'center');

end