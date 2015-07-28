4:55 PM 7/20/2015

Functions, scripts and data are mixed in the same folder. File names are exactly as found in original code.
I have most of the time made minor changes such as parameters. Sometimes I added graphs and changed code 
that did not run.

Scripts are named with S_*

================
TOPICS
================

(1) Estimation & data work

    MeanCovarianceRobustEstimators
    FancyBetaCalculation
    MissingData
    OutliersDetection
    FactorsOnDemand          -- general procedure portfolio factor algo
    PCAstatarbExample        -- Link btw OU and PCA to find co-integrated time series
    DimentionReductionBasics -- Few examples of using PCA and regression to reduce the estimation dim  
	
(2) General tools

    FullyFlexProbability -- Regime dependent probabilities into hitorical/simulated views
                            (generalization of EWMA)
    EntropyPooling       -- General methodology to input different views in historical PDF
                            (kind of generalization of Fully flex prob)
    PanicCopula          -- General way of building copula from empirical/simulated PDF
    PCAstatarbExample	 -- Applied PCA to cointegration with "trading" example

(3) Stress tests
    
    StressTestWithEntropyPooling 
    LiquidityRiskMarketImpact

(4) Allocation

    MeanVarianceOptimizationBasics   -- Basic Examples of MV
    MeanVarianceOptimizationBayesian -- Bayesian MV example
    BlackLittermanModelBasics
    EntropyPooling                   -- Generalization of Black Litterman
    DynamicAllocationStrategies      -- Basic mostly illustration of multi horizon opt and also stop-loss type of idea
    DynamicEntropyPooling            -- Multi horizon. Generalizes EntropyPooling and DynamicAllocationStrategies
    ManagingDiversification          -- Algo spreads risk uniformly given contraints (generalization of 1/n)

(5) Market Impact
    
    LiquidityRiskMarketImpact        -- Also talks about MI theory and opt execution

=======================
Bootcamp approx daily
=======================

(1) Mostly first 3 days

    MeanCovarianceRobustEstimators -- day 3
    FancyBetaCalculation           -- day 3
    MissingData                    -- day 3
    OutliersDetection              -- day 3
    FactorsOnDemand                -- day 2 & 5       
    PCAstatarbExample              -- day 2        
    DimentionReductionBasics       -- day 2 & 3
    PanicCopula                    -- day 2 & 4 (example of copula)

(2) Mostly last 3 days

    StressTestWithEntropyPooling     -- day 4
    FullyFlexProbability             -- day 4
    ManagingDiversification          -- day 4
    MeanVarianceOptimizationBasics   -- day 5
    MeanVarianceOptimizationBayesian -- day 5&6
    BlackLittermanModelBasics        -- day 6
    EntropyPooling                   -- day 6
    DynamicAllocationStrategies      -- day 6
    DynamicEntropyPooling            -- day 6 but w/o detail
    LiquidityRiskMarketImpact        -- day 6 but just mentioned

