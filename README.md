# Beta

This set of code was developed for perturbative Beta function analysis. Specifically, it employs the pade approximation with the hopes of finding roots close to that of the existing five loop order Beta function. Note that this analysis is confined to a model with SU(3) symmetry identical in most aspects to Quantum Chromodynamics (QCD). Only the number of flavors are allowed to change. Roots of this beta function could have implications for alternate gauge theories that include novel explanations for the Higgs mass. At this point, explanation of the physics is a bit out of the scope of this readme, however, I plan to include my presentations and Latex document on the subject for a little context. In this project I have attempted to include all my script with the caveat that some parts have not been completed or were scrapped early because they were not needed. I developed all this code using a BBEdit IDE and command line without any version control or external solutions.

QCD_2.py:
  The first function to run tells most of the story. It was developed procedurally and showcases the graph of coupling vs. energy as a function of coupling. To date modern physicists have only analytically derived the Beta function to loop order 5. All contemporary Beta approximations are located in the graph to the right of the screen. An adjustable slide allows the user to freely change the number of flavors. For reference the standard model includes six flavors. Two more sections were created to the left that include a visual and numerical representation of the complex zeros. The only reason for this was to determine where the zeros converge on the number line as the flavor number changes. We can easily see interestingly that they would singularly meet between models with numbers of flavors equal to 12 and 13.
  
QCD12-13.py:
  Since an integer number of flavors does not provide a concrete zero of the beta function, a natural follow up question is where exactly does the zero lie? This window was created identically to QCD12-13.py in order to find it. Physically a fractional flavor does not make sense of course, however, maybe there are brave souls who would ponder such a meaning. It is interesting to see that with modern knowledge the flavor number that is closest to zero is 12.89 and closer to 13. I think most would suspect some sort of symmetry here so more questions for further research remain. will perturbative or non perturbative methods move that root closer to 12 or 13? Anyway, this is thought provoking but programmatically identical to the above functionality.
  
Beta_Head.py:
  Enter OOP. I started by simply coding the mathematical underpinnings. As I added more and more, I encapsulated reusable peices of code into functions. When I started copying all or parts of functions to be used in seperate programs I identified the need for objects. This script is the Beta object that I built for the next stage of research.
  
Beta_Graph.py:
  This is the culmination of my research and it took me a while to get here but it was well worth it. Typically perturbative methods are done by power series and broken down into a polynomial. Mathemeticians are always trying new things and one of those hardly used power-series-like approximations is the Pade Approximation. Long story short it is two power series in a ration format and it gets fairly complicated when you try to analyze approximation schemes with higher denominator order than numerator order and vice versa all while ranging the number of flavors. (Thank God I didn't try to set them up against all five loop orders... could you imagine SU(N) symmetry?... what a combinitorics nightmare...(Challenge accepted?)) Pade approximations have a sort of unstable assymptotic behaviour in certain regions which helped me narrow down some of the more interesting approximates. Also, the Pade approximations are bounded by the original 5 loop power series. That was a long walk for a ham sandwich. Pretty much I chose three renditions of the Pade approximations to lay against the 5 loop series. Interestingly the pade approximation favors an upswing at 12 flavors to create a zero, however it is unclear if the model breaks down due to the asymptotic nature of the approximation, or if asymptotic freedom is forfeited here. 
  
***This next section features intermediary steps, unfinished or scrapped programs***

QCDbfunc.py:
  This program was one of the early additions of my programming the mathematical underpinnings for the beta function itself. I displayed my results with matplotlib and then used it as a sandbox to test slider functionality.
  
padeApprox.py:
  Coding the Pade Approximate was one of the most challenging parts of this whole process. It was a new mathematical concept that I had never heard of, and I was given the task of translating it into code. It doesn't seem extremely hard in retrospect, but it sent me bouncing between the IDE, my personal white board, and the Wolfram Alpha's article on the subject (https://mathworld.wolfram.com/PadeApproximant.html). I definitely had to rethink of the topic in matrix form rather than polynomial form for the creation of the algorithm which wasn't super easy. I had a lot of fun! This function was used to vary the M and L values for the taylor series of e^x which I could check on Wolfram Alpha to be sure my algorithm worked properly.
  
betaPade.py:
  Eventually of course I had to connect my Beta function work with my Pade Approximate work, and this took some care but I used this function to connect them and finally the code was copied into the appropriate programs.
  
BPgui.py:
  This function was actually part of my original idea for the "Beta Graph" portion of my research. I programmed it to be more versatile and able to take in any L and M values, coupling constant and number of flavors. I had plans to then be able to graph something similar to Beta_Graph.py with the use of a button. A few time consuming problems arose while I was working on this. First, there are interdependencies when picking these numbers, so there were faulty values in some regions. Also, programmatically I was spending a lot of unecessary time trying to link the graph and slider capability to the orignal window flavor number. There were a lot of little problems here and there that I may have been able to fix, however, bbecause of time constraints and  with coaching from my advisor, I was able to rethink about what I needed and save myself from getting in to deep.
  
 ***The Senior_Thesis.pdf is my paper on the project***
 
 

