debug notes:

rs values:
optimal for 497: 198
optimal for 498: -317

taking exp(model_output):
optimal for 497: 198 (0 loss)
optimal for 498: towards +inf (20.6 loss)

but apparently we can't have both at the same time :(
my reasonning: since increasing rs always improve the loss for 498, the bias and weight will keep pushing 

I found something:

I don't know if you can visualize it, but it has a strong asymptotical. If you're less than that asymptotical, you have negative output, and decreasing the input makes the output closer to 0. If the input is larger than the asymptotical, the output is positive, increasing it makes it go to 0 towards infinity and decreasing the input close to 0 makes the output go to +infinity.
To be simple, if the model predicts a data point to be on the wrong side of the asymptotical, where the optimum for this datapoint is on the other side, the model will never be able to find that optimum through gradient descent. 

Chat has some nice ideas:
1. Understanding the Problem
Your equation has the form:

𝑓
(
𝑥
)
=
𝐴
𝐵
+
𝑥
f(x)= 
B+x
A
​
 
Asymptote: At 
𝑥
=
−
𝐵
x=−B, the denominator 
𝐵
+
𝑥
B+x approaches zero, and the output 
𝑓
(
𝑥
)
f(x) diverges to 
+
∞
+∞ or 
−
∞
−∞, depending on the sign of 
𝐴
A.
Gradient Dynamics: When 
𝑥
x is on the wrong side of the asymptote, the gradient with respect to 
𝑥
x drives the output in the wrong direction, reinforcing the error.
If the model outputs 
𝑥
x such that:

𝑥
<
−
𝐵
x<−B: Gradient descent may reduce 
𝑥
x, pushing it further from the correct region.
𝑥
>
−
𝐵
x>−B: Gradients behave as expected and can optimize the loss.
The asymptote acts as a "barrier," and the model struggles to cross it during training.
