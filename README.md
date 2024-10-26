# Numerical Methods Console Application

## Overview

This project is a console application that demonstrates various numerical methods for solving different types of mathematical problems. Developed collaboratively by our group, this project provides an opportunity to understand and implement important numerical algorithms. This console based project can be used to solve systems of linear equations, non-linear equations, differential equations, and also perform matrix inversion. 

The aim of this project is to get an idea about how the computers or calculators operate based on various type of methods. It should be noted that some of the methods are better on some certain criteria. So, depending on this criteria, some equations may converge quickly than others.
Let's take a quick review on this application. 

## Team Members

- Member 1: Tawhidul Hasan
- Member 2: Oitijya Islam Auvro
- Member 3: Mamunar Rahman

## Project Structure

The application is organized into several sections, each dedicated to solving a different category of problems:

### 1. Solution of Linear Equations:

*Jacobi Iterative Method*: 

This iterative method solves a system of linear equations by repeatedly updating estimates of the solution until it converges within a maximum tolerance(0.00001 in this case). The precondition is that the elements of the matrix should be diagonally dominant. If that's not the case, then swap_rows function makes sure that the swapping of rows create a diagonally dominant matrix. 
Setup: 
Begin with an initial guess for the solution values, often starting with zeros.
Iterative Process: 
Update each variable one by one, using only the values from the previous iteration. This helps refine each solution by working independently. 
Convergence Check: 
Compare the latest values to the previous ones. If they’re close enough (within a small threshold), stop. If not, repeat until they are.

*Gauss-Seidel Iterative Method*: 

This is an advanced method over the Jacobi method, this method uses the latest values available to improve convergence. This also gets converged only if the matrix of the elements is diagonally dominant. This method takes less time than Jacobi method because of using the latest value to improve convergence. 
Setup: 
Initialize with an initial guess (usually 0), as with the Jacobi method.
In-place Update: 
For each variable, update using the most recent values (including changes from the current iteration as we need to take care of the next possible values). This accelerates convergence.
Convergence Check: 
Continue iterating until the values stabilize within a small difference.

*Gauss Elimination*: 

This method transforms the matrix into an upper triangular form, making it easier to solve using back-substitution. This is helpful in solving the equations more quickly and efficiently. By transforming into upper triangular form, we move to Gauss-Jordan elimination for converting the elements of the matrix into value 0 ( except main diagonal) . Then for each equation we transform the matrix into row echelon form and finally we get our answer for each variable. But the concerning part is that the matrix should be square matrix. That means if the unknown variable is numbered 'n' then the matrix should be in 'n×n'. Based on this condition, we can solve as many equations as possible. 
Transform to Upper Triangular: 
Rearrange the system by eliminating terms below the main diagonal, which simplifies solving.
Back Substitution: 
Once transformed, start from the last equation and work upwards, solving each equation to find the variables. In this way, the back calculation is done and we get the desires values.

*Gauss-Jordan Elimination*: 

An extension of Gauss Elimination, this method reduces the matrix further to row echelon form for a more straightforward solution. This method is also applicable for square matrix.
Reduce to Row Echelon Form: 
Start with Gauss elimination, create an upper triangular matrix. Then continue beyond Gauss elimination, eliminating terms both below and above each pivot. This creates a matrix that has non- zero values in the main diagonal but zero in other position. 
Extract Solution: 
Once reduced, each row represents a single variable directly, making it easy to identify each solution.

*LU Factorization*: 

This method decomposes a matrix into a lower and upper triangular matrix, which can simplify solving equations, especially for multiple right-hand sides.
Decompose the Matrix: 
Split the main matrix into two simpler ones: one with all elements below the main diagonal and one with elements above.
Solve in Two Stages: 
First solve for an intermediate set of values using the lower matrix, then use the upper matrix to solve for the variables.

Each method is implemented to handle systems with a minimum of 5 equations, as required.


### 2. Solution of Non-linear Equations


 *Bisection Method*: 

This method repeatedly divides an interval in half to locate a root.
Identify Interval: 
Start with an interval where the function changes sign, meaning there’s a root somewhere between.
Halve the Interval: 
Check the midpoint. If it’s close enough to zero, it’s the root; if not, adjust the interval and repeat.


 *False Position Method*: 

This method uses linear interpolation between interval endpoints to approximate the root.
Set Up Interval: 
Similar to bisection, use an interval where the function changes sign.
Approximate Root:
 Use a line between interval endpoints to estimate the root location. Adjust the interval based on this estimate until reaching the root.


*Secant Method*: 

This method uses secant lines drawn between successive approximations to approximate roots without requiring a derivative.
Initial Points: 
Start with two points near the expected root.
Iterate with Secant Approximation: 
Use a line between the points to approximate the root, updating with each iteration.


*Newton-Raphson Method*: 

This method uses the derivative of the function to rapidly converge to a root.
Start with Initial Guess: 
Begin with an estimated root.
Refine the Estimate: 
Update the guess by factoring in how steeply the function changes at that point. Repeat until it converges close to zero.

### 3. Solution of Differential Equations

*Runge-Kutta Method*: 

Specifically, the fourth-order Runge-Kutta method is implemented here. This method provides accurate solutions for first-order ordinary differential equations by estimating intermediate values in each step.
Estimate Intermediate Steps: 
Calculate several intermediate values, using them to get a highly accurate average rate of change.
Update Solution: 
Use this average to get the next value in the sequence. Continue across the entire range of values.

### 4. Matrix Inversion

This section includes a method to calculate the inverse of a matrix using matrix operations. Techniques such as elimination or LU decomposition are used to achieve this accurately.

## Getting Started
### Prerequisites
- *Programming Language*: C++

### Installation
1. Clone the project repository:
   ```bash
   git clone https://github.com/AuvroIslam/Numerical-Lab-Project.git
   cd Numerical-Lab-Project
