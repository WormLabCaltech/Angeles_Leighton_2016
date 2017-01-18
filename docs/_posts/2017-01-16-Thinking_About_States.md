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

Not quite. If you have enough data, one of the things you can do is throw darts
at your PCA, so to speak. Imagine taking many datasets as you vary only one
parameter (cell cycle, for example) and plotting them on your PCA plot. Your
final PCA plot should hopefully show some predictable pattern as the cells
vary in their cell cycle. In this way, you can begin to understand what each
axis means, although the conclusions have to be really very soft and cautious.
Nevertheless, the point remains: The PCA coordinates can in principle be
decoded and understood in a biologically relevant manner. So what? Well, arguably
the next component might be to treat the PCA coordinates as quantitative values
of some complex---yet understandable---trait.

Once this trait has been defined, then we can draw correlations between the
gene expression level and this trait of interest. The genes that are most correlated
in one direction or the other may be hypothesized to be causally related.
For an example of what I have very shoddily described here, see *this* paper.

Ok, well, but this isn't particularly easy to do, and it's not trivial to
understand. It's also got serious problem---how do you know that your trait is
really what you think it is? What does it mean to change by one unit along one
of these axes? What are the units of the axes?

Another way to try to understand this data is to give up on trying to predict the
function of every single gene in this scenario, and instead try to say something
about the relationship between the variables that we played with in this study.
Ah! Fortunately, this problem is very old, and it has been tackled very aggressively.


to be continued....
