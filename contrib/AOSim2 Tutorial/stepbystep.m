% A simple script set up as a step by step guide to AOSim2
% 
% Author: Alexander Rodack
% Date: 2/16/2015
% Software provided by Johanan L. Codona
clear all; clc; close all;

%% A Note About Using This Guide
% While the run button can be pressed to execute this entire tutorial, it
% is recommended that you read through it, and copy/paste the code into the
% command window to execute it.  This will allow you to follow along with
% the comments and execute the code in a step by step fashion.  If you want
% to run sections of it as a whole at one time, that is also possible by
% executing the "Run Section" button with the piece you want to run as the
% active cell.  For those of you that do not know, a section of Matlab
% script that starts with a double % and has a bold title is known as a
% cell.  The cell that has a yellow background (in the default color
% scheme) is the active one.  Clicking the Run Section button will execute
% the commands only in this active cell.

%% Getting Started!
% It is important to keep the updated version of AOSim2 as the version you
% are using.  The easiest way to do this is to make use of command line git
% commands.  GUIs also exist if that floats your boat, but they can be
% somewhat of a crutch.  It is always good to learn a new skill, and using
% command line is a powerful way to do things.
% 
% This is easier on Linux than on Windows, and I have no Apple computers,
% but I assume it would be closer to Linux than Windows.  This will provide
% some simple instructions based on Linux, but all the git commands will
% work no matter where you use them as long as you have git installed on
% your computer and are using a Git Shell. I do know from experience that
% you can install a git shell in Windows that works as the command line 
% tool directly from the github website.
% 
% Open the command line tool (Terminal in Linux).
% If you don't have git, in Linux type: sudo apt-get install git
% If you have another OS, go to github.com and look for more instructions.
% Or just google Git Shell. I am sure you will find something official.
% 
% If/once you have git, you are ready to make a clone of the AOSim2/OPTI528
% branch.
% 
% To do this, simply type:
% 	git clone -b OPTI528 https://github.com/jlcodona/AOSim2.git
% 
% This will create a directory with the name AOSim2. If you want to specify
% the name of the directory and where it is located, add a path to the end
% of the command above.  
% An example of this is: git clone -b OPTI528 https://github.com/jlcodona/AOSim2.git /home/alex/Desktop/AOSim_git/
% In this example, the OPTI528 branch is cloned onto the Desktop into a folder named AOSim_git
% For you, replace the /home/alex/Desktop/AOSim_git/ with <path/to/folder/name> of your desiring.
% 
% After you type in that command, press enter, and you will some lines
% appearing in Terminal.  There shouldn't be any errors, and the last line
% will read "Checking connectivity...done."
% 
% BOOM! You now have AOSim2 on your computer in the folder and location you specified.
% 
% To examine this a bit, use the cd command to switch into the directory
% that was just created.  Now type: git branch
% You should see "* OPTI528"  This means everything worked correctly.
% Now type: git status
% You will see: On branch OPTI528. Your branch is up-to-date with 'origin/OPTI528'.
% This means you have the most recent version of the repo.
% 
% The AOSim2 tutorial is in the contrib folder, and then in AOSim2 Tutorial, and is called stepbystep.m
% The example Dr. Codona went over in class on 2/17/2015 is in the examples folder.
% 
% NOTE: Make sure to periodically update your git folder from the github repository.
% Just type:
% git pull 
% 
% and you're up to date!
%
% Once you add a file (a script to run something using AOSim2) to this
% directory, you will need to add it to the repo.  Do this using the
% command: git add <path/to/file>  example: git add stepbystep.m
% Type in git status after doing this, and you should see in green
% "new file:   stepbystep.m"
% This change now needs to be committed.
% Type: git commit. A window will open for you to type in a comment about
% what has changed in the files being commited. Save that file, and the
% commit will go through.  using git commit -a will automatically commit
% all modifications (you will still need to do the changes file.  If you
% have your own git repo, you can sync to it, and push these changes so it
% is always up to date.  To do that, look up how to sync to a remote, and
% use the command: git push.  If you have code that fixes a bug, or is
% something John finds interesting, he will likely have you push it to the
% OPTI528 branch for everyone to be able to use/enjoy.
% 
% Now that you have the software, we will examine how to use it.  There
% will be more on git later.


%% Object Oriented Matlab
% For those of you who are not familiar with Object Oriented Programming
% (OOP), this will provide a little insight.  Mathworks has some good
% tutorial videos (I know from experience) that will walk you through the
% basics of creating a class structure, and using it. I highly recommend
% them if you have never seen OOP before.
%
% Basically, in OOP, you create objects that are of user defined class
% structures.  These are essentially data structures that know how to use
% the data within them to do something.  
%
% This comes to fruition in the classes by means of properties and methods.
% The properties are the data items that store whatever you want (scalars,
% vectors, matrices, cells, strings, doubles, logicals, EVERYTHING that
% is/was/can be done in Matlab).  The methods are the functions that tell
% the class what to do with that stored data.
%
% If you look at a class in AOSim2, or in the tutorial videos, you will see
% that it is broken up into sections.  It will start with the name of the
% class and the superclass structure it is going to inherit its initial
% properties and methods from.  This can be handle, matlab.mixin.Copyable,
% or the name of any user defined class (you will see this a lot in AOSim2)
% Next is a list of properties. There are protected and public properties.
% If you want to learn more, look into the tutorials.  Following the
% properties list are the methods.  The first method is always called the
% Constructor. This method tells the class how to create an object of its
% class type.  This is the place to start when looking at a new class so
% you can understand how to create an object to use in your scripting.  All
% the methods that follow (sometimes called utilities) add functionality to
% the class.  This is the part that makes the data structure know what to
% do with the data stored within it.  Finally, there are overloaded methods
% (methods that have the same name as a method in a superclass to change
% them to fit the current class) and static methods. Again, if you want to
% know more, do the tutorials. They are helpful. Really.
%
% One more note, which you would know if you watched the tutorials or are
% already in the know about OOP. You call a property or method by using a
% dot operator.  This is why you will see a lot of g.name or g.coords.
% Basically, if the property or method is known to the class, you can
% access it by using the dot.

%% Let's Jump In! Creating the (AO)Grid
% Alright, now that some of the logistical stuff is out of the way, lets do
% some stuff.
%
% AOGrid is located in the @AOGrid folder in the AOSim2 section of the
% repository.  The reaon for all the '@' symbols in front of the folder
% names is that it allows for other functions to be included in the folder
% and act as methods in the class, but also maintain the ability to be
% called separate from the class.
%
% Now that you have AOGrid open, take a moment to read the comments at the
% beginning.  These are some notes from John that are important to using
% all of AOSim2 correctly.  Also, John has some great comments, so read
% them.  Some are helpful in understanding what is happening, others are
% good for a laugh or two.
%
% Now look over the properties.  This are things that will be a part of
% every class in AOSim2 (every other class has AOGrid as a superclass).
% Some are fairly straightforward, others might be somewhat harder to
% understand what they are for.  Knowing the properties of AOGrid is less
% important than for some of the other classes, but you should get a basic
% understanding of them.  We will learn more about them as we run into the
% need to use them.  Now take a look at the Constructor method.  You can
% see here how to create an AOGrid object, and the possible input arguments
% that are accepted by the code.  Finally, take some time to look at all
% the other methods that are included.  There again is no real reason to
% explain them all here, and a lot of them are used when others are called.
% I will go into detail about some of them when the need arises.
%
% So let's make a simple AOGrid object named Grid:
Grid = AOGrid(64);

% This creates an AOGrid object that contains a 64x64 array.  The default
% values for data properties seen in the Constructor method are set to the
% object as well. 

% Let's explore a couple of the methods that can quickly become important
% coords and COORDS
% Lets run them and see what we get:
[x,y] = Grid.coords;
[X,Y] = Grid.COORDS;

% Take a look at the results
% You will see that x,y are vectors, and X,Y are meshgrid-like matrices.
% These commands map a real life coordinate system (in meters) to the
% pixels in the AOGrid based on the spacing property (defaults to 0.04)
% and how many pixels you gave to the array when creating the object (64
% in this example script).
% This is very convenient for things like plotting, and for understanding
% what a simulation might be like in physical units.

% grid
% This is another important method. This is used a lot in other classes
% that come down the line, but is useful to look at here.  Calling
% Grid.grid will print the current array stored in the object Grid. If you 
% have a matrix already created somewhere that is the right size, you can 
% set the property grid_ (the actual array is stored here) to that matrix.  
% You can also set grid_ by using the constant method.

% The following code will print what is stored in Grid.grid_ under 3
% different cases:
figure(1)
subplot(1,3,1)
imagesc(Grid.grid)
caxis([-1,1]);
title('Default Array');
subplot(1,3,2)
matrix_A = magic(64);
Grid.grid(matrix_A);
imagesc(Grid.grid)
title(sprintf('Array Set Using grid Method\n'));
subplot(1,3,3)
Grid.constant(1);
imagesc(Grid.grid)
caxis([-1,1]);
title('Array Set Using constant Method');

% Printed in the figure 1 are now the default array values (zeros), the
% array when it is set by the grid method to a matrix that is already
% known, and when it is set by the contant method.  Notice that the array
% is overwritten each time, because AOGrid can store only 1 array at a
% time.
input('Press a Key to Continue');

% Finally, take note of the overloaded operators section.  This is defining
% how objects with the class type of AOGrid react to mathematical operator
% commands (+, *, and -).  These operators become very important once we
% get to later classes (AOField especially).

% Spend some time playing around with the Grid and learn what the class is
% capable of doing.  I found that the easiest way to not only get a feel
% for using OOP, but also to begin to understand the power of AOSim2.  Call
% some methods, give inputs to them, and see what happens.  I will give you
% a few of my favorites to get you started:

% Grabbing Values from a matrix at certain coordinate points and making
% them a vector (You will see the power of this when we get to Deformable
% Mirros)
OPL = magic(64);
Grid.grid(OPL);
pistonvec = Grid.interpGrid(x(1,:),y(1,:));

% Create a circle, use Grid to fourier transform it, and then plot the
% complex result without a Cdata error.
R = sqrt(X.^2 + Y.^2);
cyl = R<=0.7;
Grid.grid(cyl);
fgrid = Grid.fft(64);
Grid.grid(fgrid);
figure(2)
subplot(1,3,1)
Grid.plotC(1);
title('OOP Designation');

% This example brings me to an important point about AOSim2.  Whenever you
% use a script that makes use of the .fft command (either from a different
% method, or straight as done in this example), the fftgrid_ property is
% set and cached.  This means you have to clear the property if you want to
% do a different FFT. It is written this way to save time from having to do
% the same FFT over and over, so if the property is not empty, the program
% assumes you want to use the cached copy and avoid unnecssary computation.
% Clearing the property is simple.  Use the touch method.
Grid.touch;

% Another note.  You don't have to use the OOP designation to call methods.
% They are all written such that the first input is the class object, which
% means you can get the same effect by using standard Matlab notation
clear fgrid;
grid(Grid,cyl);
fgrid = fft(Grid,64);
grid(Grid,fgrid);
subplot(1,3,2)
plotC(Grid,1);
title('Standard Matlab Function Designation');
touch(Grid);

% You can see in figure 2 that calling the methods either way nets the
% exact same result.  Pretty cool, right?

% Now we should take a look at the property FFTSize.  This is an important
% tool for taking Fourier Transforms in AOSim2.  This is a grid property
% that defaults to the size of the array stored in grid_.  Setting it to a
% number that is larger than the size of grid_ will result in a higher
% resolution FFT.  Lets take a look.

clear fgrid;
Grid.grid(cyl);
Grid.FFTSize = [1024,1024];
fgrid = Grid.fft; %note not giving an input to fft uses the size stored in FFTSize
subplot(1,3,3)
plotCAmpl(fgrid,0.5); %this is a function in the utils folder that functions the same way as plotC
title('Higher Resolution FFT');
Grid.touch;
% Always remember to use touch on the object after doing a FFT.

% Keep looking through the methods in AOGrid and have some fun! The more
% you know about this class, and how to use the methods for the object, the
% better off you will be as we add complexity and move forward. 

input('Press a Key to Continue');
%% Step 2: Make a Pupil for Your System
% This step can vary in complexity widely based on what you need to
% simulate.  If you just need a single, circular pupil, this is very
% simple.  That being said, one of the great things about AOSim2 is that it
% allows its users to create pupils of almost arbitrary complexity, from
% the GMT pupil to a segemented hexagonal pupil with piston/tip/tilt
% control over each segment (like the JWST or an IrisAO pupil). You just
% have to know how to use the software!
% 
% Lets start with the AOSegment class.  As you can see, the superclass of
% this is AOGrid.  Just to reiterate, this means that an AOSegment object
% inherits all the properties and methods from the AOGrid class, and adds
% in some of its own.  Take a moment to look through the class before
% moving on to the next part.
%
% Now that you have an idea of what is included in AOSegment, you might
% have some questions about it.  So lets go ahead and make one and examine
% how it works to try to straighten some thing out.
% To start, we have to define a few things:

% Pupil Vector
% PUPIL_DEF1 = [0 0 8.4 1 0.05 0 0 0 0 0];
% Let me take a second to explain this vector. Each input holds a piece of
% information that will be used when making the Pupil.
% PUPIL(1) = x0, the x-coordinate of the center of the pupil
% PUPIL(2) = y0, the y-coordinate of the center of the pupil
% PUPIL(3) = D, the diameter of the pupil
% PUPIL(4) = Transmission Type

% "Transmission" type Options:
% 0: circular hole
% 1: mirror
% -1: wedge
% -2: spider vanes (3: width, 6: nvanes, 7: starting theta)

% PUPIL(5) = smoothing or "antialiasing"
% PUPIL(6) = apodization type

% Supported Apodizations: apod_type = pupils(6)
% 0: Cosine (arg 5 holds the width)
% 1: Sonine (nu is arg 7)  (beware of definitions!)
% 2: Elliptical: arg7 is Dy.
% 3: Angel-Cheng-Woolf
% 4: Angel2
% 5: Spergel2D: (arg 5 is gaussian equivalent diameter).
% 6: Woolf slab: (arg 5 is the height)
% 7: Specified circular mask. (set shadingMask and shadingRadii
%     mask) arg 7 is an overall complex phasor.)
% 8: ellipse: (use AOAperture/addEllipse)

% If you want to learn more about all these options, go to the @AOSegment
% folder and open make1.m. make1 is the method that actually constructs the
% pupil based on the information you provide in the definition vector. It
% contains all the code that is used to draw the pieces you want.  It would
% take too long to go over every option here, so I will leave that for you
% to figure out if you happen to need to use one of them.  Instead, we will
% construct some basic pupils to help you understand the basics of what to
% do.

% The above vector PUPIL_DEF1 creates a circular pupil with nothing special
% about it.  It has a diameter 8.4 with a smoothing of 0.05.  We will see the
% effect of smoothing once we look at the pupil.  

% Now lets construct a Pupil.
Seg1 = AOSegment;
Seg1.name = 'Circular Pupil'; %this names the segment
D = 8.4; %The diameter is always in units of meters
SPACING = D/100;  %this sets the number of pixels across the segment.  It
                  %is important when it comes to sampling and numerical
                  %artifacts.
SMOOTHING = SPACING/3;


% Set up the Segment and make it
PUPIL_DEF1 = [0 0 D 1 SMOOTHING 0 0 0 0 0];
Seg1.spacing(SPACING);
Seg1.pupils = PUPIL_DEF1;
Seg1.make;

% Seg1 is now the pupil.  Lets take a look at it
figure(3)
subplot(1,3,1)
Seg1.show;
colormap(gray);
% show is a method used to plot the array stored in an object. It is
% originally a method in AOGrid, but you will notice it is overloaded in
% most of the classes to make sure it plots the most interesting aspect.

% Now lets make some other pupils of increasing complexity.
Seg2 = AOSegment;
Seg2.name = 'Obstructed Circular Pupil';
D = 8.4;
Secondary = D/5; %diameter of a secondary mirror
SPACING = D/100;
SMOOTHING = SPACING/3;
PUPIL_DEF2 = [0 0 D 1 SMOOTHING 0 0 0 0 0
              0 0 Secondary 0 SMOOTHING/2 0 0 0 0 0];
Seg2.spacing(SPACING);
Seg2.pupils = PUPIL_DEF2;
Seg2.make;
subplot(1,3,2)
Seg2.show;

Seg3 = AOSegment;
Seg3.name = 'Obstructed Circular Pupil with Spiders';
D = 8.4;
Secondary = D/5; %diameter of a secondary mirror
Spider = 0.1;
numSpiders = 4;
spider_length = D/1.9;
angle_start = 45; %degrees
SPACING = D/100;
SMOOTHING = SPACING/3;
PUPIL_DEF3 = [0 0 D 1 SMOOTHING 0 0 0 0 0
              0 0 Secondary 0 SMOOTHING/2 0 0 0 0 0
              0 0 Spider   -2 SMOOTHING numSpiders angle_start spider_length 0 0];
Seg3.spacing(SPACING);
Seg3.pupils = PUPIL_DEF3;
Seg3.make;
subplot(1,3,3)
Seg3.show;
input('Press a Key to Continue');
% As you can see, to add new pieces to the pupil, another line is included
% in the pupil definition vector.  Like I said, you can make just about any
% pupil you can imagine.

% On that note, lets show you that is actually possible by creating a
% picture of a pupil in or outside of Matlab and loading it into the grid_
% array of an AOSegment Object.
img = imread('doubleSlit.png');
% AOGrid will only use the first plane if it has more dims.
figure(3);
subplot(1,3,1)
Seg1.name = 'Loaded in Pupil';
Seg1.grid(img).show;
% Note that if you say Seg1.make it will overwrite your manual grid with a
% rendered one based on the pupils definitions.  It will not do this unless
% you tell it to.
% Now we can redefine D based on the loaded in pupil
D = max(Seg1.extent)

input('Press a Key to Continue');
% Now that we have an idea of how to define some pupil segments, lets use
% them to show the AOAperture object.  

% Lets use Seg3 to make a Binary type Telescope Pupil.  First, lets copy
% Seg3 using the deep copy method based on the superclass of AOGrid.
Seg4 = Seg3.copy;
% This method copies the information in Seg3 to a new variable Seg4.  This
% is the only way to create a new object based on an old one.  Simply using
% the command Seg4 = Seg3; just copies the pointers to Seg3.  This means if
% you change something in Seg4, you also change it in Seg3.  Utilizing the
% above method (Seg4 = Seg3.copy) creates new pointers to a new object.

% Now lets make the aperture
A = AOAperture;
A.name = 'Binary Telescope';
% Rename the segments
Seg4.name = 'Right Pupil';
Seg3.name = 'Left Pupil';
% Add each segment to the seglist of A.  The argument is which segment, and
% the offset to put it at in coordinate space (meters).
A.addSegment(Seg3,[0,-4.5]);
A.addSegment(Seg4,[0,4.5]);
% Now lets look at it.
figure(4)
A.show;

% This process can be repeated for any number of segments into a single
% aperture to combine them into a single pupil. There are a lot of methods
% in AOAperture that allow you to do things to segments individually.  This
% includes things like adding piston/tip/tilt, allowing something like an
% IrisAO Deformable Mirror to be constructed (**NOTE** if you do this, the
% tip/tilt will be relative to the pupil, not each segment, which is a bug
% that should be worked out)

input('Press a Key to Continue');

%% Now that we have a Pupil, Lets Shine Some Light on it
% AOSim2 uses the class AOField to construct objects to model electric
% fields.  You will again notice this class is a subclass of AOGrid. This
% class also includes a section of properties called Static Constants.
% These are all textbook values of wavelengths for different bandwidths. The
% cool part about these is that they can be accessed without the need of
% making an AOField Object.  Lets demonstrate that:
lambda = AOField.VBAND
% Look in the Command Window and see that lambda has the value 0.5556e-6.
% This has NOT made a field object and stored it in the name lambda. lambda
% is just a scalar. This is also separate from the AOField property lambda.
% The property is the wavelength of the electric field that is the AOField
% object. This might sound confusing, so lets make a field and see that
% they are different.

F = AOField(A);
% The input argument here is A so that the field is constructed with the
% same parameters as the Aperture we have already made (things like
% spacing and array size). This is the best way to define a field to make
% sure everything stays numerically compatible. That being said, it is
% usually a good idea to make the field object's array size larger than
% your aperture by some amount.
F.lambda = AOField.HeNe_Laser;
F.lambda
F.resize(512);

% Now lets set FFTSize to a higher value to get a better resolution when we
% do an FFT (this will be done in the mkPSF method).
F.FFTSize = 1024;

% Next we need to define the Diffraction Angle in arcseconds for later
THld = F.lambda/D * 206265;

% Now that we have all of this defined, we can do some stuff with our field
% First, lets use the planewave method to set F to be a planewave.  After
% that (in the same line) we are going to have it be incident on the
% Aperture. To do this, simply multiply them together. You might be tempted
% to dot multiply here (.*) but do not do this. You will see in AOField
% that the mtimes method has again been overloaded.  This is possibly the
% most important redefinition of this in AOSim2. Just by using *, we can
% apply phase to the field via AOScreens, AOAtmos, and DMs.
F.planewave * A;
% When we do this, Matlab will wonder why (as can be seen here with the
% underline under the * that says "'*' produces a value that might be
% unused".  You can ignore this.  The grid_ property of F has been updated.
% Now we can look at what F looks like.
figure(5)
F.show

% This isn't very interesting yet because the planewave has zero phase. You
% can see however that more space is seen around the aperture corresponding
% to the fact that we made the array size of F larger than the array size
% of A.  For those of you who are used to Matlab getting very upset when
% trying to do operations with matrices of mismatched dimmensions, this
% should make you feel like all the problems in the world aren't going to
% happen when using AOSim2. Well, that's not true. They will still happen.
% But at least in this case (and all cases of "multiplying" field objects
% by other AOSim2 class objects), this is not a problem you will see.

% This is a good stepping stone into making AOPhaseScreens to make this
% picture a little more interesting, as well as investigating the mkPSF and
% plotPSF methods in the AOField class.













