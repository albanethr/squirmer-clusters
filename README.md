# squirmer-clusters
Code for "Hydrodynamic interactions between squirmers near walls" (article under review)
All files run with Matlab. 


I/ A single squirmer above a rigid boundary
• fig3.m: This function is used to obtain the phase diagram for a single squirmer above a boundary in the far field in Fig.3. It finds the equilibrium states if they exist by splitting between puller and pushers, and finding the real roots of the relevant polynoms (eq 1.9) that are located above the wall (h > 1) if they exist.

II/ Two-dimensional symmetric encounter of two squirmers in the far field
• fig6.m: This function describes the dynamics of two approaching squirmers starting from D at their single squirmer equilibrium. It uses the auxiliary function xspeed.m , yspeed.m , and rotationrate.m and returns the plot (b,c,f,g) of figure 6 for the chosen range of β and initial separation D.
• xspeed.m: Returns the horizontal translation speed (x-direction) for two squirmers of strength β at height h above a no-slip wall, and distance D from one another.
• yspeed.m: Vertical translation speed (y-direction).
• rotationrate.m: Rotation rate for the two interracting squirmers.

III/ Near-field orientation stability of a circle of squirmers
• fig9.m Stability of circles of pullers above one or two walls uses the function isstable.m to assess the stability of a given configuration, and uses matlab function fminbnd to find the minimal value of |β| that is stable.
• fig10.m Stability of circles of pullers above one or two walls.
• isstable.m This is the auxiliary function used by fig9.m and fig10.m to determine whether a given circle, with
parameters β, ε, h is stable when bounded by nwalls walls.
• circlenearfield.m This function displays the dynamics of a circle of squirmers in the near field. It opens a figure that shows the time-evolution of the orientation of each swimmer as an arrow which goes from black to yellow with increasing time. The simulation stops if the circle is unstable and displays the time at which it broke apart
• filledcircle.m This function displays the dynamics of a circle of squirmers with a central squirmer, as used in Sec3.(g).
• partialcircle.m This function displays the dynamics of a partial circle of squirmers, as used in Sec3.(g).


