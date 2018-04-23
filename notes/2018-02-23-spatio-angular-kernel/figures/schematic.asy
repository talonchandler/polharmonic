usepackage("amsmath");
settings.outformat = "pdf";
settings.render=16;// Draw other lines
settings.embed=true;
settings.prc=true;
import three;
import graph3;
import solids;
import grid3;
defaultpen(fontsize(20pt));
currentprojection = orthographic(5,5,1.5);

texpreamble("
\usepackage{amssymb}
\usepackage{amsmath}
");

size(6cm,0);

// Draw the lower lens.
real r = 1;
real alpha = pi/6;
real f(real x) { return sqrt(r^2 - x^2); }
path s = graph(f, 0, r*sin(alpha), operator..);
path3 b3 = path3(s, plane=YZplane);
surface solidsurface = surface(b3, c=O, axis=Z);
draw(solidsurface, grey+opacity(0.25));
//label(L=Label("$\Omega$"), 0.9*r*Z);

// Draw the upper.
real r = 1;
real alpha = pi/6;
real f(real x) { return 4-sqrt(r^2 - x^2); }
path s = graph(f, 0, r*sin(alpha), operator..);
path3 b3 = path3(s, plane=YZplane);
surface solidsurface = surface(b3, c=O, axis=Z);
draw(solidsurface, grey+opacity(0.25));
//label(L=Label("$\Omega$"), 0.9*r*Z);

// Draw the polarizer
pen thinblack = black+1;
path planeoutline = box((-0.5, -0.5), (0.5, 0.5));
draw(shift(2*r*Z)*surface(planeoutline), surfacepen=grey+opacity(0.25));

// Draw the detector
pen thinblack = black+1;
path planeoutline = box((-0.5, -0.5), (0.5, 0.5));
draw(shift(4*r*Z)*surface(planeoutline), surfacepen=grey+opacity(0.25));


//Draw Axes
real ax_len = 0.5;
draw(L=Label("$\mathbf{\hat{x}}$", position=Relative(1.1), align=W), O--ax_len*X,thinblack, Arrow3(emissive(black))); 
draw(L=Label("$\mathbf{\hat{y}}$", position=Relative(1.1), align=E), O--ax_len*Y,thinblack, Arrow3(emissive(black))); 
draw(L=Label("$\mathbf{\hat{z}}$", position=Relative(1.0), align=N), O--ax_len*Z,thinblack, Arrow3(emissive(black)));

// Draw focal length labels with double arrow
triple A = 0.5*X-0.5*Y;
draw(L=Label("$f_o$", position=Relative(0.5), align=W), A--(A + r*Z), thinblack, Arrow3(emissive(black)));
draw(L=Label("$f_o$", position=Relative(0.5), align=W), (A + r*Z)--A, thinblack, Arrow3(emissive(black)));
draw(L=Label("$f_o$", position=Relative(0.5), align=W), (A+r*Z)--(A + 2*r*Z), thinblack, Arrow3(emissive(black)));
draw(L=Label("$f_o$", position=Relative(0.5), align=W), (A + 2*r*Z)--(A+r*Z), thinblack, Arrow3(emissive(black)));
draw(L=Label("$f_t$", position=Relative(0.5), align=W), (A + 3*r*Z)--(A+2*r*Z), thinblack, Arrow3(emissive(black)));
draw(L=Label("$f_t$", position=Relative(0.5), align=W), (A + 2*r*Z)--(A+3*r*Z), thinblack, Arrow3(emissive(black)));
draw(L=Label("$f_t$", position=Relative(0.5), align=W), (A + 4*r*Z)--(A+3*r*Z), thinblack, Arrow3(emissive(black)));
draw(L=Label("$f_t$", position=Relative(0.5), align=W), (A + 3*r*Z)--(A+4*r*Z), thinblack, Arrow3(emissive(black)));
draw(O--A, dashed);
draw(r*Z--(A+r*Z), dashed);
draw(3*r*Z--(A+3*r*Z), dashed);

// // Draw alpha
// real r = 1;
// real q = alpha; //theta
// real f = -pi/4; //phi
// triple B = r*expi(q,f);
// draw(O--B, dashed);
// draw(L=Label("$\alpha$", (0.5,0,1.5)), arc(O,0.15*length(B)*Z,0.5*B), thinblack);

// Choose vector
real rv = 0.5;
real q = 0.75*pi; //theta
real f = 0*pi; //phi
triple Av = rv*expi(q,f);

// Draw other lines
draw(L=Label("$\mathbf{\hat{s}}_o = \cos\varphi_o\sin\vartheta_o\hat{\mathbf{x}} + \sin\varphi_o\sin\vartheta_o\hat{\mathbf{y}} + \cos\vartheta_o\hat{\mathbf{z}}$", position=Relative(1.1), align=SE), O--Av, thinblack, Arrow3(emissive(black))); 
//draw(2*Z--(2*Z+0.5*Y), thinblack, Arrow3);
//label(L=Label("$\mathbf{\hat{p}_{\text{exc}}}$"), position=1.8*Z+0.5*Y);


// draw(L=Label("$\mathbf{\hat{s}}_o = \cos\varphi_o\sin\vartheta_o\hat{\mathbf{x}} + \sin\varphi_o\sin\vartheta_o\hat{\mathbf{y}} + \cos\vartheta_o\hat{\mathbf{z}}$", position=Relative(1.1), align=SE), O--Av, thinblack, Arrow3(emissive(black))); 
//draw(2*Z--(2*Z+0.5*Y), thinblack, Arrow3);
//label(L=Label("$\mathbf{\hat{p}_{\text{exc}}}$"), position=1.8*Z+0.5*Y);

// Plane labels
label(L=Label("$\mathbf{r}_o =x_o\hat{\mathbf{x}} + y_o\hat{\mathbf{y}}$"), align=E, -A + O);
label(L=Label("$\mathbf{r}_b = r_b\cos\phi_b\hat{\mathbf{x}} + r_b\sin\phi_b\hat{\mathbf{y}}$"), align=E, -A + 2*r*Z);
label(L=Label("$\mathbf{r}_d = r_b\cos\phi_b\hat{\mathbf{x}} + r_b\sin\phi_b\hat{\mathbf{y}}$"), align=E, -A + 4*r*Z);

// Refractive index labels
label(L=Label("$n_o$"), 0.5*A + 0.5*r*Z);
label(L=Label("$n_t$"), 0.5*A + 1.5*r*Z);
label(L=Label("$n_t$"), 0.5*A + 2.5*r*Z);
label(L=Label("$n_t$"), 0.5*A + 3.5*r*Z);

shipout(scale(4.0)*currentpicture.fit());
