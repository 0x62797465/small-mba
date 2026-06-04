# What is this
Small MBA proof of concept. This works under a prime finite field, unlike traditional MBA, which means it's also hard to apply. This can turn an expression like `a+b` into `(6448*(a)+52680*(b)+3702*((9096*(a)+12039*(b)+15981*(a^b)+47552*(a&b)+15592*(a|b)+42443*(a+b)+33192*(a-b)+47824+50782*(~a)+52862*(~b)))+38903*(a&b)+31499*(a|b)+46192*((63784*(a)+47329*(b)+22494*(a^b)+18922*((45427*(a)+22383*(b)+37549*(a^b)+2925*(a&b)+58868*(a|b)+2356*(a+b)+13423*(a-b)+60881+26581*(~a)+42212*(~b)))+39455*(a|b)+37983*(a+b)+24883*(a-b)+63312+57556*(~a)+56856*(~b)))+39721*(a-b)+26530+62040*(~a)+28830*(~b))`.
# How
You essentially create a matrix with a bunch of random values, and in the neighboring columns, you have operations applied to those values (under the finite field). By sampling the nullspace of this matrix, you get a combination of operations that leads to zero under the finite field (for all numbers, usually). All you have to do is add whatever intended operation you want to that null vector, and you've essentially created a hard-to-understand vector of operations that maps to what you originally intended. 
# Implementation
Bad. Most of the code was stuff I knew in theory, but had no idea how to implement algorithmically. While it does work, it likely has some conceptual issues that fit with each other to produce a good final result.
# To-do
A lot. Speed can be dramatically improved, obfuscation does not reorder operations, the not operation acts weird, etc. 
