# Educational_Simulations
These are simulations meant to depict the working principles behind Magnetic Particle Imaging as teaching tools. In the process these scripts make a number of simplifying assumptions (e.g. ideal, perfectly homogeneous coils, relaxation-free Langevin mangetization, etc.) but the essence of the techniques remain in tact. So far, there are three simulation functions: 
* GeneralMPIPhysics - This function is inspired by Fig 1 in Gleich and Weizenecker's 2005 paper. This function animates that figure in and (hopefully) provides slightly more detail to allow more a more intuitive  understanding of the physics

* X_Space2D_Demo - This function takes in an object (in the form of an image) and demonstrates an X-Space MPI simulation and reconstruction of that image

* SystemMatrix2D_Demo - This function takes in an object (in the form of an image) and demonstrates a System Matrix MPI simulation and reconstruction of that image
