GLIST = [1.01:.1:2 2.2:.2:5]
for ng=1:length(GLIST)
    GAIN=GLIST(ng)
    COUNTS_ = 0;
    for n=1:3000
        TSERIES = randn([1,1000]);
        BASIS=100+cumsum(TSERIES*3+.12);
        ETF=exp(cumsum(log((BASIS(2:end)./BASIS(1:end-1)-1)*GAIN+1)));
        [C,B]=hist(real(ETF./(BASIS(2:end)/BASIS(1))),BINS);
        COUNTS_=COUNTS_+C;
        %plot(BINS,cumsum(COUNTS_/sum(COUNTS_)));grid;ylim([0 1]);drawnow;
    end
    RETURNS(:,ng) = cumsum(COUNTS_/sum(COUNTS_));
end
