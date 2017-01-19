---
layout: page
title: Sperm-loss causes entry into the female state, regardless of how sperm is lost
permalink: /female_state/
---

How are *fog-2* and natural sperm depletion related? Formally, there are four
possibilities.

* *fog-2* and sperm depletion through time are independent of each other
(they do not affect the same genes at all)
* *fog-2* and sperm depletion have additive effects. They act on the same
phenotype but are otherwise independent of each other.
* *fog-2* and sperm depletion share an epistatic interaction, which means they
act along the same pathways. The combined effect of having a *fog-2* mutation
and being a 6d old adult is less than expected by additivity.
* *fog-2* and sperm depletion share a synthetic interaction, which means they
act along parallel pathways. The combined effect of having a *fog-2* mutation
and being a 6d old adult is more than expected by additivity.

In order to explore these possibilities, we decided to first ask whether
these two perturbations were acting on overlapping gene modules. If so, do their
effects have the same sign?

The plots below will help us answer these questions. We have plotted the
biological effect of each variable versus another. Points are colored according
to the logarithm of their q-value. Purple indicates a q-value closer to 0.1,
whereas yellow indicates a q-value $<10^{-20}$.

# *fog-2* partially phenocopies changes in expression due natural sperm depletion

Indeed, the graph below shows that not only do the *fog-2* transcriptome and the
transcriptome of 6 day old *C. elegans* act on the same modules, they also share
the same sign!

<iframe src="{{ site.baseurl }}/plots/aging_vs_genotype.html"
    style="max-width = 100%"
    sandbox="allow-same-origin allow-scripts"
    width="100%"
    height="600"
    scrolling="no"
    seamless="seamless"
    frameborder="0">
</iframe>
# Changes due to sperm depletion are not additive with changes in expression associated with mutation of *fog-2*

Next we checked whether the interaction coefficient was positively or negatively
correlated with the aging coefficient. A negative correlation indicates that,
on average, mutation of *fog-2* and aging are epistatic to each other and act
along the same molecular pathways; whereas a positive correlation indicates that
on average mutation of *fog-2* and aging have a synthetic interaction and they are
acting along parallel pathways.

We can see that the correlation is negative, which indicates that *fog-2* and
aging are acting along the same molecular pathways.

<iframe src="{{ site.baseurl }}/plots/aging_vs_interaction.html"
    style="max-width = 100%"
    sandbox="allow-same-origin allow-scripts"
    width="100%"
    height="600"
    scrolling="no"
    seamless="seamless"
    frameborder="0">
</iframe>
For completeness, we also show a plot of the genotype coefficient versus the
interaction coefficient. Having seen the plot above, we expect a negative
correlation between genotype and interaction. Indeed, that is what we see.
This evidence suggests that sperm-depletion, whether natural or induced by a
mutation, causes entry into a new *C. elegans* state.

<iframe src="{{ site.baseurl }}/plots/genotype_vs_interaction.html"
    style="max-width = 100%"
    sandbox="allow-same-origin allow-scripts"
    width="100%"
    height="600"
    scrolling="no"
    seamless="seamless"
    frameborder="0">
</iframe>
