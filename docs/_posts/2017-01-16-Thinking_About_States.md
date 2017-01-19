---
layout: post
title: Thinking About States - How to Use RNA-Seq in an Entirely New Way
---

Hi! Welcome to our site. I hope if you're here, it means that you enjoyed the
paper, and want to hear some more about our work. In this brief blog post, we
would like to tell you a little bit about our thoughts on how to use RNA-seq
to identify internal states in multicellular organisms.

RNA-seq, as you will know by now, is a fantastic tool with which to explore
genome-wide relationships. A frustration that comes with these tools, however,
is the feeling that they are purely a descriptive apparatus. Get your favorite
mutant, compare it to wild-type, select a new target and study that! And let's
ignore all the other hundreds of changes that we got from that study.

Well, descriptive is not always bad. But descriptive can be boring when there
are no alternatives. Luckily, this is changing. In the very recent past, many
labs, including us, have begun to use RNA-seq as a quantitative tool to discover
not gene targets, but rather genetics interactions (see *Aviv Regev*'s wonderful
paper on T-cells or maybe our BioRxiv *paper* on reconstructing the hypoxia
pathway). As you may have realized, though, the focus of *this* particular paper
was not to identify genetic interactions. Rather, we tried to identify a novel
state of the *C. elegans* life cycle. Somewhat amazingly, the statistical method
that we used to identify this state is exactly the same method that Aviv's team
used to reconstruct genetic interactions in T-cells!

I think there's something really interesting in that. Maybe I can explain, at the
risk of being terribly dull about numbers, why the formalism is the same.

First, though, a bit of history. Usually, when people do RNA-seq, they make
what is called a PCA plot. In few words, Principal Component Analysis takes an
enormous matrix of data and understands it to find the multidimensional lines
that your data lies on most compactly. Then, we plot the two lines that most
compactify your data. That gives you a compact description of the data in a way
that is easy to visualize. This is great, because plotting N-dimensional data is
a real problem. On the other hand, what do these coordinates mean? The new
coordinates are composite coordinates -- they are a number made by adding very
many components together in a weighted manner. Well, it seems we've just traded
one problem for another!


First, though, a bit of history. Usually, when people do RNA-seq, they make
what is called a PCA plot. In few words, Principal Component Analysis takes an
enormous matrix of data and understands it to find the multidimensional lines
that your data lies on most compactly. Then, we plot the two lines that most
compactify your data. That gives you a compact description of the data in a way
that is easy to visualize. This is great, because plotting N-dimensional data is
a real problem. On the other hand, what do these coordinates mean? The new
coordinates are composite coordinates -- they are a number made by adding very
many components together in a weighted manner. Well, it seems we've just traded
one problem for another! Indeed, the problem of how to interpret PCA coordinates
is rather complicated, but certain papers (for an example, see
[this one](http://www.sciencedirect.com/science/article/pii/S0092867416311497))
have begun to make progress in using PCA to understand how cellular
states are related to gene expression.

Well then.

Another way to try to understand this data is to give up on trying to predict the
function of every single gene in this scenario, and instead try to say something
about the relationship between the variables that we played with in this study.
Ah! Fortunately, this problem is very old, and it has been tackled very aggressively.

Whenever we are trying to assess the effect of a variable on an output, the first
thing to try if the trend looks linear (or if the variable is dichotomous) is
to fit a linear model. Simply put, a linear model works as follows:

Imagine that you are walking on a straight line. Suppose we are measuring the
effect that mutations on two genes $G_1$ and $G_2$ have on you. We find that
knocking out $G_1$ makes you walk two steps forward; whereas knocking out
$G_2$ makes you walk three steps forward. Both genes affect the same phenotype
(your movement), and therefore they are not entirely independent. How independent
are they? In genetics, there are three kinds of genetic interactions: non-interaction
(or what I call true independence), additive interaction (most people still call
this independence), and non-additive interaction. To test whether these two genes
interact, we should make a double mutant. If the genes are interacting additively,
then you should walk two steps + three steps = five steps. If the genes are not
interacting additively, we expect that you will walk less than five steps if the
genes are in the same pathway, or more than five steps if they share a synthetic
interaction.

People have been doing exactly this kind of analysis with qPCR for a very long
time. More recently, people have begun to do exactly this kind of analysis with
RNA-seq and long story short (if you are seriously interested, you should look
out for our other paper coming out soon) it works well. It's exactly like doing
qPCR, except 20,000 times.

But there's something else that is wonderful about this idea. So far, I've been
talking about interactions between **genes**. But there's nothing, absolutely
nothing in the book that says that the model couldn't mix interactions between
**genes** and the **environment** or **time** or **internal states**. And that is
what makes, in my opinion, this paper so special. This is the first time that we
have observed a developmental stage in *C. elegans* using transcriptomes, and the
way in which we observed it was exactly by using a linear model that measures
deviations from additivity.

So what do you think? Is this cool? Let us know!

---David Angeles
