%% Copyright (C) 2008 Simone Zuccher
%%
%% This file is NOT part of Octave but you can redistribute it and/or modify it
%% under the terms of the GNU General Public License as published by
%% the Free Software Foundation. For more details
%% see the GNU General Public License at
%% <http://www.gnu.org/licenses/>.
%%
%% This is distributed in the hope that it will be useful, but
%% WITHOUT ANY WARRANTY. 

%% -*- texinfo -*-
%% @deftypefn {Function File} {} gps
%% @deftypefnx {Function File} {} gps (@var{cmd})
%% @deftypefnx {Function File} {} gps (@var{ps-type},@var{name})
%% @deftypefnx {Function File} {} gps (@var{action}, @var{data}, @var{attribute}, @dots{})
%% Motivation/purpose.
%% Octave is a great environment to run computations but it might look not 
%% very flexible to people that use Gnuplot (http://gnuplot.sourceforge.net/)
%% for making highly-customized
%% figures, i.e. figures with special labels, arrows, mathematical 
%% formulae in the legend or on the axis, and so on.
%% A simple example is, from version 3.0.0 of Octave, the impossibility 
%% to make plots with bullets (or filled circles), available in previous versions
%% of Octave.
%%
%% Description.
%% 'gps' (GnuPlot Stream) 
%% produces two-dimensional and three-dimensional plots using a direct
%% gnuplot stream of commands sent to an X11 window. 
%% Therefore, all gnuplot features are usable.
%% The only drawback is that data passed to gps are written on local (hidden)
%% files. This is necessary to allow the full usage of 'replot'.
%% 
%% Many different combinations of arguments are possible.
%%
%% @noindent
%% No arguments.
%% @example
%% gps
%% @end example
%% This causes gps to close previously opened plots, i.e.
%% if no windows were open, nothing is done; if a gps window was open, then it 
%% is closed. Usage:
%%
%%
%%
%% @noindent
%% One arument. 
%% @example
%% gps (@var{cmd})
%% @end example
%% Possibilities for @var{cmd} are:
%% @itemize @bullet
%% @item
%% A single gnuplot command or multiple gnuplot commands separated by 
%% semicolon. Examples include
%% @example
%% gps("plot [-2*pi:2*pi][-1:1] sin(x); set title 'My first plot'")
%% gps("set xlabel 'x'");
%% gps("replot cos(x); set title 'sin and cosine of x'")
%% @end example
%%
%% As you can see, make sure you use correctly " and '.
%% Note that gps adds "replot" to every gnuplot command to make sure
%% that all listed commands are executed.
%%
%% @item
%% If @var{cmd} is "exit" or "e" or "quit" or "q", gps closes the gnuplot
%% stream, i.e. this is equivalent to invoking gps without arguments. 
%% All following examples close the gnuplot stream.
%% @example
%% gps
%% gps("exit")
%% gps("e")
%% gps("quit")
%% gps("q")
%% @end example
%%
%% @item
%% If @var{cmd} is "ps" or "pss" or "psbk" or "psbk3d" or "psqv",
%% a postscript figure named gps.ps is created from the current plot.
%% Command "ps" refers to a figure with the standard gnuplot settings, whereas
%% command "pss" refers to square axis (i.e. x and y axis have the same length).
%% Options with "bk" ("psbk" and "psbk3d") refer to the standard size(s)
%% used in the book by Squassina and Zuccher (2008).
%% Option "psqv" refers to quantum vortex figures.
%% Notice that "ps" and "pss" generate postscript 'cropped' figures, i.e. 
%% dvips tries to minimize the margins. 
%% Therefore, figures might have different BoundingBox.
%% When "psbk" or "psbk3d" or "psqv" are invoked, the BoundingBox is replaced
%% (so as to conform to the standard in Squassina and Zuccher (2008) in case
%% of bk).
%% Valid calls are
%% @example
%% gps("ps")     % -> standard gnuplot settings
%% gps("pss")    % -> square aspect ratio (2D plots)
%% gps("psbk")   % -> square aspect ratio scaled by 86%  (2D plots)
%% gps("psbk3d") % -> square aspect ratio scaled by 97%  (3D plots)
%% gps("psqv")   % -> square aspect ratio scaled by 90%  (2D plots)
%% @end example
%% @end itemize
%%
%% @noindent
%% Two aruments. 
%% @example
%% gps (@var{ps-type},@var{root-name})
%% @end example
%% The first argument is the type of postscript plot ("ps" or
%% "pss" or "psbk" or "psbk3d" or "psqv", as explained above), the second
%% one is the root of figure filename without extension (the extension ".ps"
%% is added automatically by gps). Valid calls are:
%% @example
%% gps("ps","fig1")     % -> fig1.ps, standard gnuplot settings, cropped
%% gps("pss","fig2")    % -> fig2.ps, square aspect ratio, cropped
%% gps("psbk","fig3")   % -> fig3.ps, square, scaled by 86%, fix size
%% gps("psbk3d","fig4") % -> fig4.ps, square, scaled by 97%, fix size
%% gps("psqv","fig5")   % -> fig5.ps, square, scaled by 90%, fix size
%% @end example
%%
%%
%% @noindent
%% More than two aruments (main purpose of gps). 
%% @example
%% gps (@var{action}, @var{data}, @var{attribute}, @var{data}, @var{attribute}, @var{data}, @var{attribute}, @dots{})
%% @end example
%% The first argument @var{action} is the plotting mode, i.e. 
%% @example
%% "plot" or "p" for 2D plots 
%% "splot" or "s" for 3D plots without surfaces
%% "splots" or "ss" for 3D plots with surfaces
%% "replot" or "rep" or "r" for any type of plot
%% @end example
%% The second argument @var{data} is the data to be plotted and, 
%% depending on @var{action}, @var{data} can be
%% @example
%% (a) @var{x}, @var{y}, @var{z} for "splots" or "ss" where @var{x} and @var{y} are vectors of 
%%     aritrary lengths (Nx for @var{x} and Ny for @var{x}) and @var{z} is a matrix 
%%     Nx by Ny which contains the values of z(i,j) at point of 
%%     coordinates x(i),y(j).
%% (b) a single matrix @var{data} N x M for all other values of @var{action},
%%     where N is the number of points to be plotted and M is the number
%%     of columns. Gnuplot uses by default the first two columns in 2D 
%%     plots and the first three columns in 3D plots.
%% @end example
%% The third argument @var{attribute} is the data attribute and 
%% includes title (which will appear in the legend), style (with lines or 
%% with points) or whatever option, exactly as in Gnuplot.
%% If more than 3 arguments are supplied, then their number must odd because
%% the 4th argument will be data again, the 5th will be their attributes, and so
%% forth. 
%% Valid calls are
%% @example
%% gps("plot",[x y],"t '1st line' w lp",\
%%      [x y.^2],"t '2nd line' w p",\
%%      [x -abs(y)],"t '3rd line' w l;\
%%      set auto;\
%%      set xlabel 'x';\
%%      set ylabel 'y'")
%% gps("replot",[x -y.^2],"t '2nd line' w lp")
%% gps("splot",[x y sin(x).*cos(y)],"t '3D line 1' w l;\
%%      set auto",\
%%      [y x x.*y],"t '3D line 2' w lp")
%% gps("splots",x,y,z,"w l t '';\
%%      set surface;\
%%      set hidden3d;\
%%      set ticslevel 0;\
%%      set xlabel 'x';\
%%      set ylabel 'y';\
%%      set label 'z= sin(sqrt(x^2+y^2))' at  0,0,1.1;\
%%      ")
%% @end example
%% As seen in these examples, the <<last>> attribute can include commands 
%% which can influence the whole plot such as
%% axis labels, axis ranges, plot title, arrows, labels, and so on.
%% This works as long as such commands are included
%% in the <<last>> attribute and are separated by a semicolon ";".
%% @end deftypefn

%% Purpose:    To avoid using Octave for making professional figures
%% Author:     Simone Zuccher
%% Created:    26 Oct 2008
%% Bug report: zuccher@sci.univr.it

function gps(varargin)
global __gp_var



% Quit if called without arguments
if(nargin==0)
   __gp_end_in__
   return
endif


% if a gnuplot stream is not open yet, open one
__gp_start_in__


if(nargin==1)
   cmd=varargin{1};
   switch (cmd)
      case {"exit", "quit", "q", "e"}
         __gp_end_in__
	 return
      case {"ps", "pss", "psbk", "psqv"}
         __gp_2ps_in__(cmd,"gps")
	 return
   endswitch
   % at this point it must be a gnuplot command...
   fprintf (__gp_var.FID, "%s;replot\n",cmd);
   fflush(__gp_var.FID);
   return
endif

% if the first argument is "ps", export plot as postscript figure 
%    and the second argument is interpreted the figure name -> name.ps
if (nargin==2)
   cmd=varargin{1};
   name=varargin{2};
   switch (cmd)
      case {"ps", "pss", "psbk", "psbk3d", "ps3d", "psqv"}
         __gp_2ps_in__(cmd,name)
	 return
   endswitch
endif

if ((nargin>2) && (mod(nargin-1,2)==0))
   pltoption=varargin{1};
   % this switch is needed only to decide about deleting previous files
   switch (pltoption)
      case {"plot", "p", "splot", "s", "splots", "ss"}
         % previous plotted files are not needed and, thus, are deleted
	 system("rm -f .__gp_plt_tmpfile*.dat");
         % set number of plotted curves to 0 
	 __gp_var.NC=0;
   endswitch
   % this is the main switch
   switch (pltoption)
      case {"plot", "p"}
	 __gp_plt_in__("plot",varargin{2:nargin})
      case {"splot", "s"}
	 __gp_plt_in__("splot",varargin{2:nargin})
      case {"splots", "ss"}
	 __gp_plt_in__("splots",varargin{2:nargin})
      case {"replot", "rep", "r"}
	 __gp_plt_in__("replot",varargin{2:nargin})
   endswitch
   return
endif

% If none of the above cases apply, then there is something wrong
disp('gps: unknown command or sintax error');
endfunction % gps



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function __gp_plt_in__(varargin)
global __gp_var

% find out the plot option and number of curves
gpoption = pltoption = varargin{1};
Ncurves = floor((nargin-1)/2);

% in case of "splots" (which does not exist in gnuplot, do differently
if(strcmp(pltoption,"splots"))
    gpoption = "splot";
    Ncurves = floor((nargin-1)/4);
endif

% final command to be sent to the stream
finalcmd=[gpoption " "];


% loop on number of curves
for i=1:Ncurves
   __gp_var.NC = __gp_var.NC+1;


   % generate a temporary file (hidden) name
   tmpname=[".__gp_plt_tmpfile" num2str(__gp_var.NC) ".dat"];
   
   % enter here only for surface plot, i.e. "splots"
   if(strcmp(pltoption,"splots")) 
      % get the data from inputs
      x = varargin{(i-1)*4+2};
      y = varargin{(i-1)*4+3};
      z = varargin{(i-1)*4+4};
      % get the lengths
      nx=length(x);
      ny=length(y);
      % open the temporary file and write data in right format for surface plot
      tmpfid=fopen(tmpname,"w");
      for ix=1:nx
	 for iy=1:ny
	    fprintf(tmpfid,"%f %f %f \n",x(ix),y(iy),z(ix,iy));
	 end
	 fputs(tmpfid,"\n");
      end
      % close data file
      fclose(tmpfid);
      % extract the attributes of the current data
      gpcmd = varargin{(i-1)*4+5};
      
      
   % enter here for regular gnuplot plots
   else 
      data = varargin{i*2};
      % if only one point is plotted, it is made a vector of 2 equal points
      if(size(data)(1)==1)
         data = [data;data];
      endif
      % save data to the file
      save("-ascii",tmpname,"data");
      % extract the attributes of the current data
      gpcmd = varargin{i*2+1};
   endif

   % add this plot to the final command
   if(i>1)
      finalcmd = [finalcmd ','];
   endif
   finalcmd = [finalcmd ' ''' tmpname ''' ' gpcmd];
%   sprintf("%s '%s' %s;",finalcmd,tmpname,gpcmd);
endfor
% wait 1/10 of a second to make sure files are written...
pause(.1)
% (re)plot data and add replot to ensure interpretation of additional commands
fprintf (__gp_var.FID, "%s; replot \n",finalcmd);
fflush(__gp_var.FID);
endfunction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function __gp_start_in__
global __gp_var
% if __gp_var does not exists create an empty one
if (isempty(__gp_var))
   __gp_var.FID=0;
   __gp_var.ON=false;
   __gp_var.NC=0;
endif

% If no gnuplot stream is open, open one and set the global variables 
if ((__gp_var.FID == 0) && (__gp_var.ON == false) && (__gp_var.NC == 0))
   __gp_var.FID=popen("gnuplot -noraise","w");
   __gp_var.ON=true;
   __gp_var.NC=0;
   fputs(__gp_var.FID,"set term x11\n");
   fputs(__gp_var.FID,"set yrange restore; plot 1./0. t '' \n");
%else
%   disp("La variabile esiste e non ho aperto niente\n")
endif
endfunction % __gp_start_in__
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function __gp_end_in__
global __gp_var
% If __gp_var exists, is not empty and is not initialized to zero, close stream
if ((exist("__gp_var") == 1) && (!isempty(__gp_var)) && ...
       (__gp_var.ON) && (__gp_var.FID > 0))
%   disp('quitting...');
   pclose(__gp_var.FID);
   system("rm -f .__gp_plt_tmpfile*.dat");
   system("rm -f __gp_2ps__tmpxxxfig*");
   __gp_var.FID=0;
   __gp_var.ON=false;
   __gp_var.NC=0;
   clear __gp_var
endif
endfunction % __gp_end_in__
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function __gp_2ps_in__(type,name)
global __gp_var


% Save current plot settings
fputs(__gp_var.FID,"save '.__gp_plt_tmpfileCURRENTPLOT'\n");
fflush(__gp_var.FID);
pause(.2);
% say why we are waiting
fdisp (stdout, ['generating figure ' name '.ps...']);
fflush(stdout);


% defaul option is regular gnuplot size
sizeoption = "noratio  1,1 ";
BBoption = "18 546 398 796";
fontoption = "12";
switch (type)
      case {"ps"}
         sizeoption = "noratio  1,1 "; 
	 BBoption = "26 547 370 790";
      case {"pss"}
         sizeoption = "square"; 
	 BBoption = "66 535 349 806";
      case {"psbk"}
         sizeoption = "square 0.86,0.86";
%         sizeoption = "square 0.9,0.9";
         BBoption = "51 567 295 809";
      case {"psqv"}
         sizeoption = "square 0.9,0.9";
         BBoption = "51 567 310 809";
      case {"psbk3d"}
         fontoption = "8";
         sizeoption = "square 0.97, 0.97";
%         BBoption = "70 550 298 743";
         BBoption = "70 560 298 753";
	 fputs(__gp_var.FID,"set format x '%g \\vspace{-3mm}\\hspace{1mm}'\n");
         fputs(__gp_var.FID,"set format y '%g \\vspace{0mm}\\hspace{-5mm}'\n");
%         fputs(__gp_var.FID,"set format z '%g \\vspace{1mm}\\hspace{-4mm}'\n");
         fflush(__gp_var.FID);
      case {"ps3d"}
         fontoption = "8";
%         sizeoption = "square 0.97, 0.97";
         BBoption = "60 570 350 775";
	 fputs(__gp_var.FID,"set format x '%g \\vspace{-3mm}\\hspace{1mm}'\n");
         fputs(__gp_var.FID,"set format y '%g \\vspace{0mm}\\hspace{-5mm}'\n");
         fflush(__gp_var.FID);
endswitch

%for i=1:2
% Export from gnuplot to pslatex
% Replace bounding box only for psbk and psbk3d
switch (type)
   case {"psqv"}
      fputs(__gp_var.FID,["set term pslatex color " fontoption "\n"]);
   otherwise
      fputs(__gp_var.FID,["set term pslatex monochrome " fontoption "\n"]);
endswitch
fputs(__gp_var.FID,"set output '__gp_2ps__tmpxxxfig.tex'\n");
fputs(__gp_var.FID,["set size " sizeoption "\n"]);
%fputs(__gp_var.FID,"set key spacing 1.7\n");
fputs(__gp_var.FID,"replot\n");
fflush(__gp_var.FID);

% Reset things as before
fputs(__gp_var.FID,"set term x11\n");
pause(.1);
fputs(__gp_var.FID,"load '.__gp_plt_tmpfileCURRENTPLOT'\n");
fputs(__gp_var.FID,"replot\n");
fflush(__gp_var.FID);
%end

% Remove temporary file
system("rm -f .__gp_plt_tmpfileCURRENTPLOT");

% Create a temporary LaTeX file to include the figure generated above
system("echo ' \\\\documentclass[12pt,a4paper]{article}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\setlength{\\\\topmargin}{-40 mm}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\setlength{\\\\oddsidemargin}{-25 mm}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\setlength{\\\\textwidth}{200 mm}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\setlength{\\\\textheight}{200 mm}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\setlength{\\\\footskip}{0 mm}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\usepackage{graphics}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\usepackage{amsmath}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\pagestyle{empty}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\begin{document}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\input{__gp_2ps__tmpxxxfig.tex}' >> __gp_2ps__tmpxxxfiglatexfile.tex");
system("echo ' \\\\end{document}' >> __gp_2ps__tmpxxxfiglatexfile.tex");

% Run LaTeX without showing anything on the screen
system("latex __gp_2ps__tmpxxxfiglatexfile.tex > /dev/null");

% Make the postscript file in silent mode
system("dvips -E -o __gp_2ps__tmpxxxfiglatexfile.ps __gp_2ps__tmpxxxfiglatexfile -q");

% Replace bounding box only for psbk and psbk3d
switch (type)
   case {"ps", "psbk", "psbk3d", "ps3d", "psqv"}
      system(["sed s/'BoundingBox: '/'BoundingBox: " BBoption " % '/ \
        < __gp_2ps__tmpxxxfiglatexfile.ps > __gp_2ps__tmpxxxfiglatexfilenew.ps"]);
      system("cp __gp_2ps__tmpxxxfiglatexfilenew.ps __gp_2ps__tmpxxxfiglatexfile.ps");
endswitch

% Change the title of the figure to the correct one and finish
systemcm = ["sed s/'__gp_2ps__tmpxxxfiglatexfile.dvi'/'" name ".dvi'/ < \
        __gp_2ps__tmpxxxfiglatexfile.ps > " name ".ps"];
system(systemcm);

% Clear up temporary files
system("rm -f __gp_2ps__tmpxxxfiglatexfile* __gp_2ps__tmpxxxfig.tex");

% Done
endfunction % __gp_2ps_in__
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
