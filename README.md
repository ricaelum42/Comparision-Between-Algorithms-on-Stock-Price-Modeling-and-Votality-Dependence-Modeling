# Comparision-Between-Algorithms-on-Stock-Price-Modeling-and-Votality-Dependence-Modeling
Winter 2018 STAT440 Final Project

## Summary
The stock of a corporation is constituted of the equity stock of its owners. A single share of the stock represents fractional ownership of the corporation in proportion to the total number of shares. Given the history price, we would like to forecast the stock prices. In time series, SARIMA is popular choice. However, if we apply SARIMA into real data, the prediction is not always as good as expected. In this project, we would like to use a diﬀerent method, decompose the stock price into three pieces, conquer each one and unit all three. 

In addition, for usual GARCH models, the marginal distributions are estimated. However, usually in real case, stocks from same industry may have their prices depend on each other. Therefore, in this project, we tried to model dependence of stock prices using copula approach. Our analysis shows that k-nearest neighbour is our best model in ﬁtting one single stock, and GARCH dependence model works best when there are more stocks in the same industry.
