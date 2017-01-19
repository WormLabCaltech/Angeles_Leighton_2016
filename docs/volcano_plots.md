---
layout: page
title: Volcano Plots
permalink: /volcano/
---

Volcano plots are a convenient way of visualizing RNA-seq data. The
x-axis is a measurement of biological effect of a variable (time, genotype
or the interaction between them). The y-axis is a measurement of the probability
that a deviation from the mean as extreme or more than the observed value
happens given that the null hypothesis is true. In practice, we are often
interested in the genes that have the lowest q-value and the largest effect
size.


# Aging Volcano

<iframe src="{{ site.baseurl }}/plots/aging_volcano.html"
    style="max-width = 100%"
    sandbox="allow-same-origin allow-scripts"
    width="100%"
    height="600"
    scrolling="no"
    seamless="seamless"
    frameborder="0">
</iframe>


# Genotype Volcano

<iframe src="{{ site.baseurl }}/plots/genotype_volcano.html"
    style="max-width = 100%"
    sandbox="allow-same-origin allow-scripts"
    width="100%"
    height="600"
    scrolling="no"
    seamless="seamless"
    frameborder="0">
</iframe>


# Interaction Volcano

<iframe src="{{ site.baseurl }}/plots/interaction_volcano.html"
    style="max-width = 100%"
    sandbox="allow-same-origin allow-scripts"
    width="100%"
    height="600"
    scrolling="no"
    seamless="seamless"
    frameborder="0">
</iframe>
