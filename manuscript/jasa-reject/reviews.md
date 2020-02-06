## Associate Editor Comments to the Author:

This is a well-written paper with an interesting application. However, both reviewers are concerned with the lack of discussion of and comparison to the existing literature. This is a serious shortcoming that must be addressed.


## Reviewer: 1 Comments to the Author

This paper proposes the use of a Markov-switching state space model to analyze musical data. The paper is well-written, and the problem is interesting. The author(s) derive the model from many specific concerns due to the domain knowledge. However, the model does not seem to be novel. It looks very similar to an existing work (Gu and Raphael, 2012), which also use a Markov-switching model to analyze musical data. The author(s) uses this model to justify the validity of their model due to similarity, but does not discuss any fundamental differences. Except simple validation or description of the fitting, the analysis of the data does not provide any additional/special insights that we gain from the proposed technique. Although I believe the technique would be a useful tool for analyzing musical data, I have concern with the novelty of the paper.

1. Is the heterogeneity in variance modeled well by the model?
2. In Algorithm 1, how sensitive is the algorithm to the choice of B? Could the author(s) provide a simple numerical study?
3. In Sec. 2.4, does it matter if we fix \sigma_{acc} and \sigma_{stress} to other values instead of 1? Will the conclusions in the real data analysis change dramatically? Is there any reason to pick 1, not other numbers?
4. Could the author(s) comment and compare with Gu and Raphael (2012)?
5. As noted by the author(s), the jumps make many smoothing techniques non-applicable/poorly performing. Could the author(s) compare with piecewise polynomial regression techniques and/or other techniques with change-point modeling?
6. In the abstract, the author(s) mention functional data analysis technique is used. But I couldn’t find any. Could the author(s) be more specific?

## Reviewer: 2 Comments to the Author

This paper presents an application of switching state space models for music performance analysis, focusing on a particular musical example (a mazurka by F. Chopin). Overall, the paper is relevant within the "Applications & Case Studies" section of JASA, and the area of musical performance analysis using statistical methodologies is an active area of research. The manuscript is very well-written and clear; model choices are always clearly justified, and the authors make a considerable effort both defending choices in the methodology as well as acknowledging limitations of the employed methods and models. Some aspects of placing this work within the context of related literature are missing, and the subsequent analysis of this particular composition is convincing, although is missing a proper evaluation (which admittedly would be hard to perform). In terms of novelty, the application of switching state space models to inferring performance attributes can be deemed as sufficiently novel. 

Detailed comments can be found below:

1. Section 1: the rather long quotes defining classical music are perhaps redundant and in any case are not a focus of the paper. A more concise definition of classical music can be directly provided by the authors, or relevant literature providing definitions can be cited. On this topic, I disagree with the authors' statement that "What separates classical music from other types of music is that the music itself is written down but performed millions of times..." - clearly there are many other music styles which involve music notation which is reproduced multiple times in various interpretations (e.g. pop/rock lead sheets).

2. Section 1.1: this section includes citations for work related to music generation and music information retrieval, which is not directly relevant to this work. There is a large body of work that focuses more specifically on topics related to this paper (i.e. statistical/mathematical methods for music performance analysis), see for example the body of work created by Profs. Elaine Chew and Gerhard Widmer (in particular the former, who has used the Mazurka dataset for dynamics-related music performance analysis).

3. Equation 1: define F.

4. Section 2.2: one aspect not really mentioned in the paper, which however feels important in the context of the music composition/performances being analyzed, is rubato. Not all music performances contain this large amount of accelerations and declarations as this particular piece which originates from a specific period where such performance practices where commonplace.

5. Section 2.2: one assumption that is not explicitly discussed in the paper is whether it's reasonable to assume that tempo decisions can be broken down into 4 discrete states, or whether there might be more nuance between those states. It would be interesting for example to observe the posterior probability of the switching states and to see whether there are multiple states being active concurrently (and what this means with respect to performance interpretation).

6. Table 3: to what extent do these priors depend on this particular Mazurka, and can these be applied to other classical music compositions?

7. Section 3.2: At this point, one would wish that there was any performance ground truth where the model's inferred states could be compared against. For example, is the lack of acceleration and stress states in Hatto's performance a due to a limitation of the model, or is this indeed an aspect of Hatto's performance? I understand that obtaining performance annotations can be difficult and also could be viewed as subjective, however some additional evaluation on the model's performance seems necessary.

8. Section 3.3: Instead of the variable-specific treatment with respect to distances, was it not possible to normalize all variables before computing the distance matrix?


