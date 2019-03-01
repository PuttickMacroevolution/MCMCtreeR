<!DOCTYPE html>

<html xmlns="http://www.w3.org/1999/xhtml">

<head>

<meta charset="utf-8" />
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<meta name="generator" content="pandoc" />


<meta name="author" content="Mark Puttick" />
<meta name="author" content="marknputtick@gmail.com" />
<meta name="author" content="University of Bath" />

<meta name="date" content="2019-01-17" />


<h1 class="title toc-ignore">MCMCTreeR</h1>
<h4 class="author"><em>Mark Puttick</em></h4>
<h4 class="author"><em><a href="mailto:marknputtick@gmail.com">marknputtick@gmail.com</a></em></h4>
<h4 class="author"><em>University of Bath</em></h4>
<h4 class="date"><em>17 January 2019</em></h4>

</div>

<ul>
<li><strong>A guide to setting-up MCMCtree analyses </strong>is shown below and a PDF is available <a href="https://github.com/PuttickMacroevolution/MCMCTreeR/blob/master/vignettes/MCMCtree.pdf">here</a></li>
<li><strong>A guide to plotting mcmc trees and showing uncertainty </strong>like the ones show below can be found as a PDF <a href="https://github.com/PuttickMacroevolution/MCMCTreeR/blob/master/vignettes/MCMCtree_plot.pdf">here</a></li>
</ul>

<div class="figure" style="text-align: center">
<img src="https://github.com/PuttickMacroevolution/MCMCTreeR/blob/master/vignettes/MCMCtree_plot_fig/plot1.png", alt="MCMCTreeR plotting" width="576" />
</div>

<p>This is a guide for the R program <a href="https://github.com/PuttickMacroevolution">MCMCTreeR</a>.</p>
<p>MCMCTreeR contains functions to set up analyses in the <a href="http://abacus.gene.ucl.ac.uk/software/paml.html">MCMCtree</a> program. The functions here help users choose the best parameters to reflect age information for prior age distributions, visualise time priors, and produce output files ready to be used in MCMCtree. A separate <a href="https://github.com/PuttickMacroevolution/MCMCTreeR/blob/master/vignettes/MCMCtree_plot.pdf">vignette</a> is available to explain the plotting options for timescaled trees in MCMCTreeR.</p>
<p><a href="http://abacus.gene.ucl.ac.uk/software/paml.html">MCMCtree</a> is a Bayesian program in the software PAML to estimate divergence times on fixed topologies using molecular data developed by Ziheng Yang. The program requires various inputs from the user: a phylogeny, molecular data, and selected model parameters. This guide does not include details about which time priors are most appropriate for the data, etc, so please see the <a href="http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf">MCMCtree manual</a>.</p>
<div id="installation" class="section level1">
<h1>Installation</h1>
<pre class="r"><code>## check for devtools and install if necessary
if (!any(rownames(installed.packages()) == &quot;devtools&quot;)) install.packages(&quot;devtools&quot;)
library(devtools)
install_github(&quot;PuttickMacroevolution/MCMCTreeR&quot;, quiet = TRUE)</code></pre>
<p>The examples here use a phylogeny of apes and associated age information. These data are a phylogeny of apes <code>apeTree</code>, the minimum ages for internal nodes <code>minimumTimes</code>, maximum ages for internal nodes <code>maximumTimes</code>, and the tip labels descending from each node <code>monophyleticGroups</code>. These example data can be substituted for other data.</p>
<pre class="r"><code>library(MCMCTreeR)
## Loading required package: ape
## Loading required package: sn
## Loading required package: stats4
## 
## Attaching package: 'sn'
## The following object is masked from 'package:stats':
## 
##     sd
## Loading required package: coda</code></pre>
<pre class="r"><code>data(apeData)
names(apeData)
## [1] &quot;minimumTimes&quot;       &quot;maximumTimes&quot;       &quot;monophyleticGroups&quot;
## [4] &quot;apeTree&quot;</code></pre>
To make these data more easily available, extract them into separate variables (as used in the code below)
<pre class="r"><code>minimumTimes <- apeData$minimumTimes
maximumTimes <- apeData$maximumTimes
monophyleticGroups <- apeData$monophyleticGroups</code></pre>
</div>
<div id="estimate-parameters-for-node-input-parameters" class="section level1">
<h1>Estimate parameters for node input parameters</h1>
<p>This section includes information to estimate and plot prior age distributions for node(s) used in MCMCtree divergence time estimation. The data required to do this is a phylogeny, minimum and maximum ages for the nodes with prior distributions, and taxa that descend from each node.</p>
<p>The code can be used to simultaneously estimate the parameter values that reflect the <strong>a priori</strong> time information for nodes and write files ready to be input in MCMCtree. MCMCTreeR can produce output files with the same type of distributions used to summarise <strong>a priori</strong> time information for all nodes, or separate distributions can be used to reflect uncertainty on different internal nodes.</p>
<p>By default, the functions here estimates the distribution parameters so that the distribution spans for user-supplied minimum bounds (lower age) and maximum bounds (upper age). By default, minimum ages are treated as ‘hard’ constraints and maximum ages are ‘soft’. The function then ensures that 97.5% of the distribution falls between these minimum and maximum ages. The code can estimate parameters for the Cauchy, Skew-T, Skew-normal, and Gamma distributions shown in the <a href="http://abacus.gene.ucl.ac.uk/software/pamlDOC.pdf">MCMCtree manual</a> on page 50.</p>
<div id="skew-t" class="section level2">
<h2>Skew T</h2>
<div id="estimate-scale-with-a-given-shape" class="section level3">
<h3>estimate scale with a given shape</h3>
<p>The default arguments in the <code>estimateSkewT</code> assumes the user wants to estimate the scale of the distribution with a given shape value (the default shape value is 50). The function estimates the parameters with the user-supplied minimum and maximum ages for all nodes, and the monophyletic groups that define the nodes. The output <code>skewT_results$mcmctree</code> shows the estimated parameters in the input ready for MCMCtree. Here the parameters for the Skew T distributions are the location (lower node age), scale, shape, and degrees of freedom</p>
<pre class="r"><code>skewT_results &lt;- estimateSkewT(minAge = minimumTimes, maxAge = maximumTimes, 
    monoGroups = monophyleticGroups, phy = apeTree, plot = FALSE)
## [1] &quot;warning - minProb parameter value recycled&quot;
## [1] &quot;warning - maxProb parameter value recycled&quot;
## [1] &quot;warning - estimateScale argument recycled&quot;
## [1] &quot;warning - estimateShape argument recycled&quot;
## [1] &quot;warning - estimateMode argument recycled&quot;
## [1] &quot;warning - shape parameter value recycled&quot;
## [1] &quot;warning - addMode parameter value recycled&quot;</code></pre>
<pre class="r"><code>skewT_results$mcmctree
##      
## 1 7 1
##                                                                                                                                                         
## 1 ((((human,(chimpanzee,bonobo))'ST(0.8,0.016,50,1)',gorilla)'ST(0.6,0.024,50,1)',(orangutan,sumatran)'ST(1.3,0.028,50,1)'),gibbon)'ST(1.5,0.059,50,1)';
##                
## 1 //end of file</code></pre>
<p>The function <code>plotMCMCTree</code> plots the estimated age distributions given these parameters.</p>
<pre class="r"><code>par(mfrow = c(2, 2), family = &quot;Palatino&quot;)
for (i in 1:4) plotMCMCTree(skewT_results$parameters[i, ], method = &quot;skewT&quot;, 
    title = paste0(&quot;node &quot;, i), upperTime = max(maximumTimes))</code></pre>
<div class="figure" style="text-align: center">
<img src="vignettes/MCMCTree_fig/plot1.png", align="middle", alt="Skew T distributions for all nodes" width="576" />
<p class="caption">
Skew T distributions for all nodes
</p>
</div>
<pre class="r"><code>skewT_results$mcmctree
##      
## 1 7 1
##                                                                                                                                                         
## 1 ((((human,(chimpanzee,bonobo))'ST(0.8,0.016,50,1)',gorilla)'ST(0.6,0.024,50,1)',(orangutan,sumatran)'ST(1.3,0.028,50,1)'),gibbon)'ST(1.5,0.059,50,1)';
##                
## 1 //end of file</code></pre>
<p>If the distributions are acceptable, the output can be written into a tree file ready to be input into MCMCtree using the function <code>estimateSkewT</code>. The functions will be written when the argument <code>writeMCMCtree</code> is set to <code>TRUE</code>, and the file is set using the <code>mcmcTreeName</code> argument. Additionally, a PDF file is output showing the estimated distributions if <code>plot=TRUE</code> and the file name can be specifying using the argument <code>pdfOutput</code>.</p>
<pre class="r"><code># result in tree MCMCTree format
skewT_results$mcmctree
##      
## 1 7 1
##                                                                                                                                                         
## 1 ((((human,(chimpanzee,bonobo))'ST(0.8,0.016,50,1)',gorilla)'ST(0.6,0.024,50,1)',(orangutan,sumatran)'ST(1.3,0.028,50,1)'),gibbon)'ST(1.5,0.059,50,1)';
##                
## 1 //end of file</code></pre>
<pre class="r"><code>## not run skewT_results &lt;- estimateSkewT(minAge=minimumTimes,
## maxAge=maximumTimes, monoGroups=monophyleticGroups,
## phy=apeTree, plot=FALSE, pdfOutput='skewTPlot.pdf',
## writeMCMCTree=TRUE, mcmcTreeName='skewTInput.tre')</code></pre>
<p>It is not necessary to specify the same shape value for each parameter: a different value of the shape parameter can be set for each distribution.</p>
<pre class="r"><code>## not run (remove ## to run) skewT_results &lt;-
## estimateSkewT(minAge=minimumTimes, maxAge=maximumTimes,
## monoGroups=monophyleticGroups, shape=c(9, 10, 8, 10),
## phy=apeTree, plot=TRUE, pdfOutput='skewTPlot.pdf',
## writeMCMCTree=TRUE, mcmcTreeName='skewTInput.tre')
## skewT_results$parameters</code></pre>
</div>
<div id="estimate-shape-with-a-given-scale" class="section level3">
<h3>estimate shape with a given scale</h3>
<p>The function <code>estimateSkewT</code> will take input minimum input times, and estimate the value of the shape that will produce the desired distribution with the scale parameter set to 0.05.</p>
<pre class="r"><code>skewT_results &lt;- estimateSkewT(minAge = minimumTimes[2], maxAge = maximumTimes[2], 
    monoGroups = monophyleticGroups, scale = 0.05, estimateShape = TRUE, 
    estimateScale = FALSE, phy = apeTree, plot = FALSE, writeMCMCTree = FALSE)</code></pre>
<pre class="r"><code>skewT_results$parameters
##        location scale shape df
## node_1      0.6  0.05     1  1</code></pre>
</div>
</div>
<div id="skew-normal" class="section level2">
<h2>Skew normal</h2>
<p>The <code>estimateSkewNormal</code> function estimates the value of the scale that will produce a skew normal distribution with the 97.5% cumulative probability of the distribution at the maximum age.</p>
<pre class="r"><code>skewNormal_results &lt;- estimateSkewNormal(minAge = minimumTimes, 
    maxAge = maximumTimes, monoGroups = monophyleticGroups, addMode = 0.05, 
    phy = apeTree, plot = FALSE)
## [1] &quot;warning - minProb parameter value recycled&quot;
## [1] &quot;warning - maxProb parameter value recycled&quot;
## [1] &quot;warning - estimateScale argument recycled&quot;
## [1] &quot;warning - estimateShape argument recycled&quot;
## [1] &quot;warning - estimateMode argument recycled&quot;
## [1] &quot;warning - shape parameter value recycled&quot;
## [1] &quot;warning - addMode parameter value recycled&quot;</code></pre>
<pre class="r"><code>skewNormal_results$parameters
##        location scale shape
## node_1     1.55  0.65    50
## node_2     0.65  0.25    50
## node_3     0.85  0.16    50
## node_4     1.35  0.29    50</code></pre>
<p>As only a single value is provided for each parameter by the user, the function outputs warnings to indicate these values are recycled for each node.</p>
<p>These skew normal distributions can be plotted to the screen.</p>
<pre class="r"><code>par(mfrow = c(2, 2), family = &quot;Palatino&quot;)
for (i in 1:4) plotMCMCTree(skewNormal_results$parameters[i, 
    ], method = &quot;skewNormal&quot;, title = paste0(&quot;node &quot;, i), upperTime = max(maximumTimes))</code></pre>
<div class="figure" style="text-align: center">
<img src="vignettes/MCMCTree_fig/plot2.png", align="middle", alt="Skew normal distributions for all nodes" width="576" />
<p class="caption">
Skew normal distributions for all nodes
</p>
</div>
</div>
<div id="cauchy" class="section level2">
<h2>Cauchy</h2>
<p>Here the <code>estimateCauchy</code> function is used to estimate parameters and plot the example on page 50 of PAML manual.</p>
<pre class="r"><code>example_page_50 &lt;- estimateCauchy(minAge = 1, maxAge = 4.32, 
    monoGroups = monophyleticGroups[[1]], phy = apeTree, offset = 0.5, 
    minProb = 0.025, plot = FALSE)[[1]]
plotMCMCTree(example_page_50, method = &quot;cauchy&quot;, title = paste0(&quot;node &quot;, 
    i), upperTime = max(maximumTimes))</code></pre>
<div class="figure" style="text-align: center">
<img src="vignettes/MCMCTree_fig/plot3.png", align="middle",  alt="Cauchy distributions for all nodes  (with a given scale)" width="576" />
<p class="caption">
Cauchy distributions for all nodes (with a given scale)
</p>
</div>
<div id="estimate-scale-with-a-given-shape-1" class="section level3">
<h3>estimate scale with a given shape</h3>
<p>The <code>estimateCauch</code> function will take input minimum input times, and estimate the value of the scale parameter that will produce a Cauchy distribution with the 97.5% cumulative probability of the distribution at the user-supplied maximum age.</p>
<pre class="r"><code>cauchy_results &lt;- estimateCauchy(minAge = minimumTimes, maxAge = maximumTimes, 
    monoGroups = monophyleticGroups, offset = 0.5, phy = apeTree, 
    plot = FALSE)
## [1] &quot;warning - minProb parameter value recycled&quot;
## [1] &quot;warning - maxProb parameter value recycled&quot;
## [1] &quot;warning - offset parameter value recycled&quot;
## [1] &quot;warning - scale parameter value recycled&quot;
## [1] &quot;warning - estimateScale argument recycled&quot;</code></pre>
<pre class="r"><code>cauchy_results$parameters
##         tL   p     c     pL
## node_1 1.5 0.5 0.075 1e-300
## node_2 0.6 0.5 0.008 1e-300
## node_3 0.8 0.5 0.001 1e-300
## node_4 1.3 0.5 0.016 1e-300</code></pre>
<pre class="r"><code>par(mfrow = c(2, 2), family = &quot;Times&quot;)
for (i in 1:4) plotMCMCTree(cauchy_results$parameters[i, ], method = &quot;cauchy&quot;, 
    title = paste0(&quot;node &quot;, i), upperTime = max(maximumTimes))</code></pre>
<div class="figure" style="text-align: center">
<img src="vignettes/MCMCTree_fig/plot4.png", align="middle", alt="Cauchy distributions for all nodes (with a given shape)" width="576" />
<p class="caption">
Cauchy distributions for all nodes (with a given shape)
</p>
</div>
<p>These plots indicate we may have constrained our distribution too much for the 2nd, 3rd, and 4th distribution so we can modify that to allow for a smaller offset.</p>
<pre class="r"><code>cauchy_results &lt;- estimateCauchy(minAge = minimumTimes, maxAge = maximumTimes, 
    monoGroups = monophyleticGroups, offset = c(0.5, 0.1, 0.1, 
        0.05), phy = apeTree, plot = FALSE)</code></pre>
<pre><code>## [1] &quot;warning - minProb parameter value recycled&quot;
## [1] &quot;warning - maxProb parameter value recycled&quot;
## [1] &quot;warning - scale parameter value recycled&quot;
## [1] &quot;warning - estimateScale argument recycled&quot;</code></pre>
<pre class="r"><code>cauchy_results$parameters</code></pre>
<pre><code>##         tL    p     c     pL
## node_1 1.5 0.50 0.075 1e-300
## node_2 0.6 0.10 0.035 1e-300
## node_3 0.8 0.10 0.022 1e-300
## node_4 1.3 0.05 0.040 1e-300</code></pre>
<pre class="r"><code>par(mfrow = c(2, 2), family = &quot;Times&quot;)
for (i in 1:4) plotMCMCTree(cauchy_results$parameters[i, ], method = &quot;cauchy&quot;, 
    title = paste0(&quot;node &quot;, i), upperTime = maximumTimes[i])</code></pre>
<div class="figure" style="text-align: center">
<img src="vignettes/MCMCTree_fig/plot5.png", align="middle", alt="Cauchy distributions for all nodes (with a given shape) and smaller offset" width="576" />
<p class="caption">
Cauchy distributions for all nodes (with a given shape) and smaller offset
</p>
</div>
</div>
</div>
<div id="uniform-distribution" class="section level2">
<h2>Uniform distribution</h2>
<pre class="r"><code>uniform_results &lt;- estimateBound(minAge = minimumTimes, maxAge = maximumTimes, 
    monoGroups = monophyleticGroups, phy = apeTree, plot = FALSE)
## [1] &quot;warning - minProb parameter value recycled&quot;
## [1] &quot;warning - maxProb parameter value recycled&quot;</code></pre>
<pre class="r"><code>uniform_results$parameters
##         tL  tU    pL    pU
## node_1 1.5 3.0 0.025 0.025
## node_2 0.6 1.2 0.025 0.025
## node_3 0.8 1.2 0.025 0.025
## node_4 1.3 2.0 0.025 0.025</code></pre>
<pre class="r"><code>par(mfrow = c(2, 2), family = &quot;Palatino&quot;)
for (i in 1:4) plotMCMCTree(uniform_results$parameters[i, ], 
    method = &quot;bound&quot;, title = paste0(&quot;node &quot;, i), upperTime = maximumTimes[i] + 
        1)</code></pre>
<div class="figure" style="text-align: center">
<img src="vignettes/MCMCTree_fig/plot6.png", align="middle", alt="Uniform distributions for all nodes" width="576" />
<p class="caption">
Uniform distributions for all nodes
</p>
</div>
</div>
<div id="gamma-distribution" class="section level2">
<h2>Gamma distribution</h2>
<pre class="r"><code>gamma_results &lt;- estimateGamma(minAge = minimumTimes, maxAge = maximumTimes, 
    monoGroups = monophyleticGroups, alpha = 188, beta = 2690, 
    offset = 0.1, phy = apeTree, plot = FALSE)
## [1] &quot;warning - alpha parameter value recycled&quot;
## [1] &quot;warning - beta parameter value recycled&quot;
## [1] &quot;warning - offset parameter value recycled&quot;
## [1] &quot;warning - estimateAlpha argument recycled&quot;
## [1] &quot;warning - estimateBeta argument recycled&quot;</code></pre>
<pre class="r"><code>gamma_results$parameters
##        alpha beta
## node_1  4304 2690
## node_2  1883 2690
## node_3  2421 2690
## node_4  3766 2690</code></pre>
<pre class="r"><code>par(mfrow = c(2, 2), family = &quot;Palatino&quot;)
for (i in 1:4) plotMCMCTree(gamma_results$parameters[i, ], method = &quot;gamma&quot;, 
    title = paste0(&quot;node &quot;, i), upperTime = maximumTimes[i])</code></pre>
<div class="figure" style="text-align: center">
<img src="vignettes/MCMCTree_fig/plot7.png", align="middle", alt="Gamma distributions for all nodes" width="576" />
<p class="caption">
Gamma distributions for all nodes
</p>
</div>
</div>
<div id="upper-age" class="section level2">
<h2>Upper Age</h2>
<pre class="r"><code>upper_results &lt;- estimateUpper(maxAge = maximumTimes, monoGroups = monophyleticGroups, 
    rightTail = 0.025, phy = apeTree)
## [1] &quot;warning - maxProb parameter value recycled&quot;</code></pre>
<pre class="r"><code>upper_results$parameters
##         tU    pR
## node_1 3.0 0.025
## node_2 1.2 0.025
## node_3 1.2 0.025
## node_4 2.0 0.025</code></pre>
</div>
<div id="fixed-ages" class="section level2">
<h2>Fixed ages</h2>
<pre class="r"><code>fixed_results &lt;- estimateFixed(minAge = minimumTimes[1], monoGroups = monophyleticGroups[[1]], 
    phy = apeTree)</code></pre>
<pre><code>fixed_results
## $parameters
## fixed age.nodeOne 
##               1.5 
## 
## $apePhy
## 
## Phylogenetic tree with 7 tips and 6 internal nodes.
## 
## Tip labels:
##  human, chimpanzee, bonobo, gorilla, orangutan, sumatran, ...
## Node labels:
## [1] &quot;'=1.5'&quot;
## 
## Rooted; no branch lengths.
## 
## $mcmctree
##      
## 1 7 1
##                                                                               
## 1 ((((human,(chimpanzee,bonobo)),gorilla),(orangutan,sumatran)),gibbon)'=1.5';
##                
## 1 //end of file
## 
## $nodeLabels
## [1] &quot;'=1.5'&quot;</code></pre>
</div>
<div id="different-parameters-on-different-nodes" class="section level2">
<h2>Different parameters on different nodes</h2>
<p>If we want different parameters on different nodes we can specify this by using the <code>mcmcTreePhy</code> function. Here there are different distributions applied to the internal nodes: a fixed root (node 1), skew normal (node 2), gamma (node 3), and upper distribution (node 4) to our tree. For each input we give the associated parameter values in a vector in the order of nodes. i.e - for the <code>min pro</code> on four nodes can be set as 1, 2, 4 to be 1e-8 and node 3 to be 0.025</p>
<pre class="r"><code>each.node.method &lt;- c(&quot;skewT&quot;, &quot;cauchy&quot;, &quot;gamma&quot;, &quot;upper&quot;)
output.full &lt;- mcmcTreePhy(phy = apeTree, minAge = minimumTimes, 
    maxAge = maximumTimes, monoGroups = monophyleticGroups, method = each.node.method, 
    writeMCMCTree = FALSE)
## [1] &quot;length of some parameters and nodes do not match - first parameter will be used for each node&quot;</code></pre>
<p>This can fine-tuned. For example, to estimate alpha not beta for the 3rd node</p>
<pre class="r"><code>estimate.alpha &lt;- c(FALSE, FALSE, TRUE, FALSE)
estimate.beta &lt;- c(TRUE, TRUE, FALSE, TRUE)
outputFull &lt;- mcmcTreePhy(phy = apeTree, minAges = minimumTimes, 
    maxAges = maximumTimes, monoGroups = monophyleticGroups, 
    method = each.node.method, estimateAlpha = estimate.alpha, 
    estimateBeta = estimate.beta, alphaInput = 188, betaInput = 100, 
    writeMCMCTree = FALSE)
## [1] &quot;length of some parameters and nodes do not match - first parameter will be used for each node&quot;</code></pre>
<p>Outputs from individual methods can be added to the input for subsequent node estimation. This allows for easier fine-tuning. Perhaps easier to explain with an example. Here, a skew normal is applied to the first node.</p>
<pre class="r"><code>skewNormal_results_nodeOne &lt;- estimateSkewNormal(minAge = minimumTimes[1], 
    maxAge = maximumTimes[1], monoGroups = monophyleticGroups[[1]], 
    addMode = 0.05, phy = apeTree, plot = FALSE, writeMCMCTree = FALSE)</code></pre>
<pre><code>skewNormal_results_nodeOne$apePhy
## 
## Phylogenetic tree with 7 tips and 6 internal nodes.
## 
## Tip labels:
##  human, chimpanzee, bonobo, gorilla, orangutan, sumatran, ...
## Node labels:
## [1] &quot;'SN[1.55~0.65~50]'&quot;
## 
## Rooted; no branch lengths.</code></pre>
<p>This output is then used as input in <code>estimateCauchy</code> to estimate parameters for a Cauchy distribution, which is applied to the second node.</p>
<pre class="r"><code>cauchy_results_nodeTwo &lt;- estimateCauchy(minAge = minimumTimes[2], 
    maxAge = maximumTimes[2], monoGroups = monophyleticGroups[[2]], 
    offset = 0.5, phy = skewNormal_results_nodeOne$apePhy, plot = FALSE, 
    writeMCMCTree = FALSE)</code></pre>
<pre><code>cauchy_results_nodeTwo$apePhy
## 
## Phylogenetic tree with 7 tips and 6 internal nodes.
## 
## Tip labels:
##  human, chimpanzee, bonobo, gorilla, orangutan, sumatran, ...
## Node labels:
## [1] &quot;'SN[1.55~0.65~50]'&quot;        NA                         
## [3] &quot;'L[0.6~0.5~0.008~1e-300]'&quot;
## 
## Rooted; no branch lengths.</code></pre>
<p>The third node is the set as a uniform distribution.</p>
<pre class="r"><code>uniform_results_nodeThree &lt;- estimateBound(minAge = minimumTimes[3], 
    maxAge = maximumTimes[3], monoGroups = monophyleticGroups[[3]], 
    phy = cauchy_results_nodeTwo$apePhy, plot = FALSE, writeMCMCTree = FALSE)
uniform_results_nodeThree$apePhy</code></pre>
<pre><code>## 
## Phylogenetic tree with 7 tips and 6 internal nodes.
## 
## Tip labels:
##  human, chimpanzee, bonobo, gorilla, orangutan, sumatran, ...
## Node labels:
## [1] &quot;'SN[1.55~0.65~50]'&quot;        NA                         
## [3] &quot;'L[0.6~0.5~0.008~1e-300]'&quot; &quot;'B[0.8~1.2~0.025~0.025]'&quot; 
## 
## Rooted; no branch lengths.</code></pre>
<p>The fourth is a skewT distribution, and the tree can be written to file for input into MCMCtree.</p>
<pre class="r"><code>## not run skewT_results_nodeFour &lt;-
## estimateSkewT(minAge=minimumTimes[4],
## maxAge=maximumTimes[4], monoGroups=monophyleticGroups[[4]],
## scale=0.5, phy=cauchy_results_nodeTwo$apePhy, plot=FALSE,
## writeMCMCTree = TRUE) skewT_results_nodeFour$apePhy</code></pre>
</div>
</div>
