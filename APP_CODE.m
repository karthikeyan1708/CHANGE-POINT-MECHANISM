classdef mech_analyser_exported_2 < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        Label                          matlab.ui.control.Label
        thetaButton                    matlab.ui.control.Button
        enterinitialangularvelocityEditField  matlab.ui.control.NumericEditField
        enterinitialangularvelocityEditFieldLabel  matlab.ui.control.Label
        PHY113Label                    matlab.ui.control.Label
        GUIFORCHANGEPOINTMECHANISMANALYSISLabel  matlab.ui.control.Label
        GROUP14Label                   matlab.ui.control.Label
        CloseButton                    matlab.ui.control.Button
        accelerationanalysisButton_2   matlab.ui.control.Button
        velocityanalysisButton         matlab.ui.control.Button
        TextArea_6                     matlab.ui.control.TextArea
        ANALYSISButton                 matlab.ui.control.Button
        TextArea_5                     matlab.ui.control.TextArea
        TextArea_4                     matlab.ui.control.TextArea
        TextArea_3                     matlab.ui.control.TextArea
        checktypeButton                matlab.ui.control.Button
        TextArea_2                     matlab.ui.control.TextArea
        TextArea                       matlab.ui.control.TextArea
        SIMULATEButton                 matlab.ui.control.Button
        clearButton                    matlab.ui.control.Button
        checkforchangepointButton      matlab.ui.control.Button
        enterlengthofr3EditField_2     matlab.ui.control.NumericEditField
        enterlengthofr3EditField_2Label  matlab.ui.control.Label
        enterlengthofr1EditField_2     matlab.ui.control.NumericEditField
        enterlengthofr1EditField_2Label  matlab.ui.control.Label
        enterlengthofr2EditField_2     matlab.ui.control.NumericEditField
        enterlengthofr2EditField_2Label  matlab.ui.control.Label
        enterlengthofr0EditField       matlab.ui.control.NumericEditField
        enterlengthofr0EditFieldLabel  matlab.ui.control.Label
        UIAxes6                        matlab.ui.control.UIAxes
        UIAxes5                        matlab.ui.control.UIAxes
        UIAxes4_2                      matlab.ui.control.UIAxes
        UIAxes4                        matlab.ui.control.UIAxes
        UIAxes3                        matlab.ui.control.UIAxes
        UIAxes2                        matlab.ui.control.UIAxes
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: clearButton
        function clearButtonPushed(app, event)
            app.enterlengthofr3EditField_2.Value = 0;
            app.enterlengthofr1EditField_2.Value = 0;
            app.enterlengthofr2EditField_2.Value = 0;
            app.enterlengthofr0EditField.Value = 0;
            app.enterinitialangularvelocityEditField.Value = 0;
            app.TextArea_6.Visible = 'off' ;
            app.TextArea_5.Visible = 'off' ;
            app.TextArea_4.Visible = 'off' ;
            app.TextArea_3.Visible = 'off' ;
            app.TextArea_2.Visible = 'off' ;
            app.TextArea.Visible = 'off' ;
        end

        % Button pushed function: checkforchangepointButton
        function checkforchangepointButtonPushed(app, event)
            r1 = app.enterlengthofr1EditField_2.Value;
            r2 = app.enterlengthofr2EditField_2.Value;
            r3 = app.enterlengthofr3EditField_2.Value;
            r0 = app.enterlengthofr0EditField.Value;
            
            [Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(r0,r1,r2,r3) ;
            assert (r1 + r2 + r3 > r0,"Not a valid Mechanism")
            assert (r0 + r2 + r3 > r1,"Not a valid Mechanism")
            assert (r1 + r0 + r3 > r2,"Not a valid Mechanism")
            assert (r1 + r2 + r0 > r3,"Not a valid Mechanism")
            
            
            function [Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(L0,L1,L2,L3)
                A = [L0,L1,L2,L3];
                [Lmax,indx_max] = max(A);
                [Lmin,indx_min] = min(A);
                i=1;
                k=1;
                while(i<=length(A))
                    if(i~=indx_max && i~=indx_min)
                        L(k) = A(i);
                        indx(k) = i;
                        k = k+1;
                    end
                    i=i+1;
                end
            end
            
            if(Lmax + Lmin == L(1) + L(2))
                app.TextArea_2.Visible = 'on' ;
                app.checktypeButton.Visible = 'on' ;
                app.SIMULATEButton.Visible = 'on' ;
                app.ANALYSISButton.Visible = 'on' ;
                app.TextArea.Visible = 'off' ;
                
            else
                app.TextArea.Visible = 'on' ;
                app.TextArea_2.Visible = 'off' ;
            end
        end

        % Button pushed function: checktypeButton
        function checktypeButtonPushed(app, event)
            r1 = app.enterlengthofr1EditField_2.Value;
            r2 = app.enterlengthofr2EditField_2.Value;
            r3 = app.enterlengthofr3EditField_2.Value;
            r0 = app.enterlengthofr0EditField.Value;
            
            if(r1==r2&&r3==r0&&r1==r3)
                app.TextArea_6.Visible = 'on' ;
                app.TextArea_4.Visible = 'off' ;
                app.TextArea_5.Visible = 'off' ;
                app.TextArea_3.Visible = 'off' ;
            
            elseif(r0==r1&&r3==r2&&r0~=r3)
                    app.TextArea_4.Visible = 'on' ;
                    app.TextArea_3.Visible = 'off' ;
                    app.TextArea_5.Visible = 'off' ;
                    app.TextArea_6.Visible = 'off' ;
                
            elseif(r1==r3&&r2==r0&&r1~=r2)
                    app.TextArea_5.Visible = 'on' ;
                    app.TextArea_3.Visible = 'off' ;
                    app.TextArea_6.Visible = 'off' ;
                    app.TextArea_4.Visible = 'off' ;
                    
            else
                app.TextArea_3.Visible = 'on' ;
                app.TextArea_5.Visible = 'off' ;
                app.TextArea_6.Visible = 'off' ;
                app.TextArea_4.Visible = 'off' ;
            end
        end

        % Button pushed function: SIMULATEButton
        function SIMULATEButtonPushed(app, event)
        

theta1 = sym('theta1');
theta2 = sym('theta2');
theta3 = sym('theta3');
r0 = sym('r0');
r1 = sym('r1');
r2 = sym('r2');    
r3 = sym('r3');


r1 = app.enterlengthofr1EditField_2.Value;
r2 = app.enterlengthofr2EditField_2.Value;
r3 = app.enterlengthofr3EditField_2.Value;
r0 = app.enterlengthofr0EditField.Value;

h = app.enterinitialangularvelocityEditField.Value;
w = rad2deg(-(h));

[Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(r0,r1,r2,r3)
assert (r1 + r2 + r3 > r0,"Not a valid Mechanism")
assert (r0 + r2 + r3 > r1,"Not a valid Mechanism")
assert (r1 + r0 + r3 > r2,"Not a valid Mechanism")
assert (r1 + r2 + r0 > r3,"Not a valid Mechanism")

assert (Lmax + Lmin == L(1) + L(2),"Not a change point mechanism")



if Lmin == r1 || Lmin == r0
    
th2=[];
th3=[];
theta1max = 0;
tmd = rad2deg(theta1max)
tmdi = 360;
time = (tmd - tmdi)/w
t=0:0.1:time

th1=( w*t + tmd)

%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3);


%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th1)

    theta1 = th1(i)
    F = @(x) [+r1*cosd(theta1) + r2*cosd(x(1)) - r3*cosd(x(2)) - r0;
             +r1*sind(theta1) + r2*sind(x(1)) - r3*sind(x(2))];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0);
    th2 = [th2 x(1)] ;
    th3 = [th3 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end

elseif Lmin == r3
    
th2=[];
th1=[];
theta3max = 0;
tmd = rad2deg(theta3max)
tmdi = -360;
time = -(tmd - tmdi)/w
t=0:0.1:time



th3=( w*t + tmd)



%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3);


%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th3)

    theta3 = th3(i)
    F = @(x) [+r1*cosd(x(2)) + r2*cosd(x(1)) - r3*cosd(theta3) - r0;
             +r1*sind(x(2)) + r2*sind(x(1)) - r3*sind(theta3)];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0);
    th2 = [th2 x(1)] ;
    th1 = [th1 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end
   
    
elseif Lmin == r2
    

th2=[];
th3=[];
h = app.enterinitialangularvelocityEditField.Value;
w = rad2deg(-(h));

theta1max = acos((r0^2 + r1^2 -  ((r3 + r2)^2))/(2*r1*r0))

tmd = rad2deg(theta1max)
tmdi = -tmd
time = -(tmd - tmdi)/w
t=0:0.1:time



th1=( w*t + tmd)



%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3)

%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th1)

    theta1 = th1(i)
    F = @(x) [+r1*cosd(theta1) + r2*cosd(x(1)) - r3*cosd(x(2)) - r0;
             +r1*sind(theta1) + r2*sind(x(1)) - r3*sind(x(2))];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0,options);
    th2 = [th2 x(1)] ;
    th3 = [th3 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end
th1
th2
th3
    
    
end



% Coordinates A B C D 
B=[]
C=[]

A=[0,0].*ones(size(th1))'
D=[r0,0].*ones(size(th1))'

for i=1:length(th1)
B=[ B ;[r1*cosd(th1(i)),r1*sind(th1(i))] ] ; 

C=[ C ;[(r0 +  r3*cosd(th3(i))),r3*sind(th3(i))]] ;

end
B
C



 i=1;
 maxlink = 2*max([r0 r1 r2 r3]);
 figure(1)
    while(i <= length(A))
    cla
    axis([-maxlink maxlink -maxlink maxlink])
        x1 =[A(i,1) B(i,1)];
        y1 =[A(i,2) B(i,2)];
        x2 = [C(i,1) D(i,1)];
        y2 = [C(i,2) D(i,2)];
        x3 = [B(i,1) C(i,1)];
        y3 = [B(i,2) C(i,2)];

        
      
        plot(0,0,"Color",'k',"LineWidth",2,"Marker","^","MarkerFaceColor","k")
        
        hold("on")
        
        plot(r0,0,"Color",'k',"LineWidth",2,"Marker","^","MarkerFaceColor","k")
        plot(x1,y1,"LineWidth",2);
        
        plot(x2,y2,"LineWidth",2);
        plot(x3,y3,"LineWidth",2);
        
        
        drawnow
        
        pause(0.1)
        
        i = i + 1;
    end

    
        

    



function [Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(L0,L1,L2,L3)
A = [L0,L1,L2,L3];
[Lmax,indx_max] = max(A);
[Lmin,indx_min] = min(A);
i=1;
k=1;
while(i<=length(A))
    if(i~=indx_max && i~=indx_min)
        L(k) = A(i);
        indx(k) = i;
        k = k+1;
    end
    i=i+1;
end
end
        end

        % Button pushed function: ANALYSISButton
        function ANALYSISButtonPushed(app, event)
            app.accelerationanalysisButton_2.Visible = 'on' ;
            app.velocityanalysisButton.Visible = 'on' ;
            app.thetaButton.Visible = 'on' ;
        end

        % Button pushed function: thetaButton
        function thetaButtonPushed(app, event)
            app.UIAxes4.Visible = 'off' ;
            cla(app.UIAxes4);
            app.UIAxes4_2.Visible = 'off' ;
            cla(app.UIAxes4_2);
            app.UIAxes3.Visible = 'off' ;
            cla(app.UIAxes3);
            app.UIAxes2.Visible = 'off' ;
            cla(app.UIAxes2);
            app.UIAxes5.Visible = 'on' ;
            app.UIAxes6.Visible = 'on' ;
            
theta1 = sym('theta1');
theta2 = sym('theta2');
theta3 = sym('theta3');
r0 = sym('r0');
r1 = sym('r1');
r2 = sym('r2');    
r3 = sym('r3');


r1 = app.enterlengthofr1EditField_2.Value;
r2 = app.enterlengthofr2EditField_2.Value;
r3 = app.enterlengthofr3EditField_2.Value;
r0 = app.enterlengthofr0EditField.Value;

h = app.enterinitialangularvelocityEditField.Value;
w = rad2deg(-(h));

[Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(r0,r1,r2,r3)
assert (r1 + r2 + r3 > r0,"Not a valid Mechanism")
assert (r0 + r2 + r3 > r1,"Not a valid Mechanism")
assert (r1 + r0 + r3 > r2,"Not a valid Mechanism")
assert (r1 + r2 + r0 > r3,"Not a valid Mechanism")

assert (Lmax + Lmin == L(1) + L(2),"Not a change point mechanism")



if Lmin == r1 || Lmin == r0
    
th2=[];
th3=[];
theta1max = 0;
tmd = rad2deg(theta1max)
tmdi = 360;
time = (tmd - tmdi)/w
t=0:0.1:time

th1=( w*t + tmd)

%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3);


%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th1)

    theta1 = th1(i)
    F = @(x) [+r1*cosd(theta1) + r2*cosd(x(1)) - r3*cosd(x(2)) - r0;
             +r1*sind(theta1) + r2*sind(x(1)) - r3*sind(x(2))];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0);
    th2 = [th2 x(1)] ;
    th3 = [th3 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end

elseif Lmin == r3
    
th2=[];
th1=[];
theta3max = 0;
tmd = rad2deg(theta3max)
tmdi = -360;
time = -(tmd - tmdi)/w
t=0:0.1:time



th3=( w*t + tmd)



%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3);


%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th3)

    theta3 = th3(i)
    F = @(x) [+r1*cosd(x(2)) + r2*cosd(x(1)) - r3*cosd(theta3) - r0;
             +r1*sind(x(2)) + r2*sind(x(1)) - r3*sind(theta3)];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0);
    th2 = [th2 x(1)] ;
    th1 = [th1 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end
   
    
elseif Lmin == r2
    

th2=[];
th3=[];
h = app.enterinitialangularvelocityEditField.Value;
w = rad2deg(-(h));

theta1max = acos((r0^2 + r1^2 -  ((r3 + r2)^2))/(2*r1*r0))

tmd = rad2deg(theta1max)
tmdi = -tmd
time = -(tmd - tmdi)/w
t=0:0.1:time



th1=( w*t + tmd)



%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3)

%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th1)

    theta1 = th1(i)
    F = @(x) [+r1*cosd(theta1) + r2*cosd(x(1)) - r3*cosd(x(2)) - r0;
             +r1*sind(theta1) + r2*sind(x(1)) - r3*sind(x(2))];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0,options);
    th2 = [th2 x(1)] ;
    th3 = [th3 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end
th1
th2
th3
    
    
end

%Velocity
for i = 1:length(th1)
    A = [ r2*cosd(th2(i)) ,- r3*cosd(th3(i));...
          -r2*sind(th2(i)) ,+ r3*sind(th3(i)) ];
    B = [-r1*w*cosd(th1(i));
        r1*w*sind(th1(i))];
    x = inv(A)*B;
    w2(i) = x(1,1);
    w3(i) = x(2,1);
end

%Acceleration
for i = 1:length(th1)
    A = [ r2*cosd(th2(i)) ,- r3*cosd(th3(i));...
          -r2*sind(th2(i)) ,+ r3*sind(th3(i)) ];
    
    B = [r1*w^2*sind(th1(i)) + r2*w2(i)^2*sind(th2(i)) - r3*w3(i)^2*sind(th3(i));...
        r1*w^2*cosd(th1(i)) + r2*w2(i)^2*cosd(th2(i)) - r3*w3(i)^2*cosd(th3(i))];
    x = A\B;
    a2(i) = x(1,1);
    a3(i) = x(2,1);
end

t1 = t;
t1(1) = [];
w2(1) = [];
w3(1) = [];
a2(1) = [];
a3(1) = [];


function [Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(L0,L1,L2,L3)
A = [L0,L1,L2,L3];
[Lmax,indx_max] = max(A);
[Lmin,indx_min] = min(A);
i=1;
k=1;
while(i<=length(A))
    if(i~=indx_max && i~=indx_min)
        L(k) = A(i);
        indx(k) = i;
        k = k+1;
    end
    i=i+1;
end
end
            plot(app.UIAxes5,t,th2,'LineWidth',2);
            cla(app.UIAxes4);
            cla(app.UIAxes4_2);
            cla(app.UIAxes3);
            cla(app.UIAxes2);
            plot(app.UIAxes6,t,th3,'LineWidth',2);
        end

        % Button pushed function: velocityanalysisButton
        function velocityanalysisButtonPushed(app, event)
            app.UIAxes3.Visible = 'on' ;
            app.UIAxes2.Visible = 'on' ;
            app.UIAxes4.Visible = 'off' ;
            cla(app.UIAxes4);
            app.UIAxes4_2.Visible = 'off' ;
            cla(app.UIAxes4_2);
            app.UIAxes5.Visible = 'off' ;
            cla(app.UIAxes5);
            app.UIAxes6.Visible = 'off' ;
            cla(app.UIAxes6);
            
theta1 = sym('theta1');
theta2 = sym('theta2');
theta3 = sym('theta3');
r0 = sym('r0');
r1 = sym('r1');
r2 = sym('r2');    
r3 = sym('r3');


r1 = app.enterlengthofr1EditField_2.Value;
r2 = app.enterlengthofr2EditField_2.Value;
r3 = app.enterlengthofr3EditField_2.Value;
r0 = app.enterlengthofr0EditField.Value;

h = app.enterinitialangularvelocityEditField.Value;
w = rad2deg(-(h));

[Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(r0,r1,r2,r3)
assert (r1 + r2 + r3 > r0,"Not a valid Mechanism")
assert (r0 + r2 + r3 > r1,"Not a valid Mechanism")
assert (r1 + r0 + r3 > r2,"Not a valid Mechanism")
assert (r1 + r2 + r0 > r3,"Not a valid Mechanism")

assert (Lmax + Lmin == L(1) + L(2),"Not a change point mechanism")



if Lmin == r1 || Lmin == r0
    
th2=[];
th3=[];
theta1max = 0;
tmd = rad2deg(theta1max)
tmdi = 360;
time = (tmd - tmdi)/w
t=0:0.1:time

th1=( w*t + tmd)

%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3);


%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th1)

    theta1 = th1(i)
    F = @(x) [+r1*cosd(theta1) + r2*cosd(x(1)) - r3*cosd(x(2)) - r0;
             +r1*sind(theta1) + r2*sind(x(1)) - r3*sind(x(2))];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0);
    th2 = [th2 x(1)] ;
    th3 = [th3 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end

elseif Lmin == r3
    
th2=[];
th1=[];
theta3max = 0;
tmd = rad2deg(theta3max)
tmdi = -360;
time = -(tmd - tmdi)/w
t=0:0.1:time



th3=( w*t + tmd)



%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3);


%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th3)

    theta3 = th3(i)
    F = @(x) [+r1*cosd(x(2)) + r2*cosd(x(1)) - r3*cosd(theta3) - r0;
             +r1*sind(x(2)) + r2*sind(x(1)) - r3*sind(theta3)];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0);
    th2 = [th2 x(1)] ;
    th1 = [th1 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end
   
    
elseif Lmin == r2
    

th2=[];
th3=[];
h = app.enterinitialangularvelocityEditField.Value;
w = rad2deg(-(h));

theta1max = acos((r0^2 + r1^2 -  ((r3 + r2)^2))/(2*r1*r0))

tmd = rad2deg(theta1max)
tmdi = -tmd
time = -(tmd - tmdi)/w
t=0:0.1:time



th1=( w*t + tmd)



%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3)

%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th1)

    theta1 = th1(i)
    F = @(x) [+r1*cosd(theta1) + r2*cosd(x(1)) - r3*cosd(x(2)) - r0;
             +r1*sind(theta1) + r2*sind(x(1)) - r3*sind(x(2))];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0,options);
    th2 = [th2 x(1)] ;
    th3 = [th3 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end
th1
th2
th3
    
    
end

%Velocity
for i = 1:length(th1)
    A = [ r2*cosd(th2(i)) ,- r3*cosd(th3(i));...
          -r2*sind(th2(i)) ,+ r3*sind(th3(i)) ];
    B = [-r1*w*cosd(th1(i));
        r1*w*sind(th1(i))];
    x = inv(A)*B;
    w2(i) = x(1,1);
    w3(i) = x(2,1);
end

%Acceleration
for i = 1:length(th1)
    A = [ r2*cosd(th2(i)) ,- r3*cosd(th3(i));...
          -r2*sind(th2(i)) ,+ r3*sind(th3(i)) ];
    
    B = [r1*w^2*sind(th1(i)) + r2*w2(i)^2*sind(th2(i)) - r3*w3(i)^2*sind(th3(i));...
        r1*w^2*cosd(th1(i)) + r2*w2(i)^2*cosd(th2(i)) - r3*w3(i)^2*cosd(th3(i))];
    x = A\B;
    a2(i) = x(1,1);
    a3(i) = x(2,1);
end

t1 = t;
t1(1) = [];
w2(1) = [];
w3(1) = [];
a2(1) = [];
a3(1) = [];


function [Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(L0,L1,L2,L3)
A = [L0,L1,L2,L3];
[Lmax,indx_max] = max(A);
[Lmin,indx_min] = min(A);
i=1;
k=1;
while(i<=length(A))
    if(i~=indx_max && i~=indx_min)
        L(k) = A(i);
        indx(k) = i;
        k = k+1;
    end
    i=i+1;
end
end

            plot(app.UIAxes2,t1,w2,'LineWidth',2);
            cla(app.UIAxes4);
            cla(app.UIAxes4_2);
            cla(app.UIAxes5);
            cla(app.UIAxes6);
            plot(app.UIAxes3,t1,w3,'LineWidth',2);
        end

        % Button pushed function: accelerationanalysisButton_2
        function accelerationanalysisButton_2Pushed(app, event)
            app.UIAxes3.Visible = 'off' ;
            cla(app.UIAxes3);
            app.UIAxes2.Visible = 'off' ;
            cla(app.UIAxes2);
            app.UIAxes4.Visible = 'on' ;
            app.UIAxes4_2.Visible = 'on' ;
            app.UIAxes5.Visible = 'off' ;
            cla(app.UIAxes5);
            app.UIAxes6.Visible = 'off' ;
            cla(app.UIAxes6);
            
theta1 = sym('theta1');
theta2 = sym('theta2');
theta3 = sym('theta3');
r0 = sym('r0');
r1 = sym('r1');
r2 = sym('r2');    
r3 = sym('r3');


r1 = app.enterlengthofr1EditField_2.Value;
r2 = app.enterlengthofr2EditField_2.Value;
r3 = app.enterlengthofr3EditField_2.Value;
r0 = app.enterlengthofr0EditField.Value;

h = app.enterinitialangularvelocityEditField.Value;
w = rad2deg(-(h));

[Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(r0,r1,r2,r3)
assert (r1 + r2 + r3 > r0,"Not a valid Mechanism")
assert (r0 + r2 + r3 > r1,"Not a valid Mechanism")
assert (r1 + r0 + r3 > r2,"Not a valid Mechanism")
assert (r1 + r2 + r0 > r3,"Not a valid Mechanism")

assert (Lmax + Lmin == L(1) + L(2),"Not a change point mechanism")



if Lmin == r1 || Lmin == r0
    
th2=[];
th3=[];
theta1max = 0;
tmd = rad2deg(theta1max)
tmdi = 360;
time = (tmd - tmdi)/w
t=0:0.1:time

th1=( w*t + tmd)

%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3);


%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th1)

    theta1 = th1(i)
    F = @(x) [+r1*cosd(theta1) + r2*cosd(x(1)) - r3*cosd(x(2)) - r0;
             +r1*sind(theta1) + r2*sind(x(1)) - r3*sind(x(2))];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0);
    th2 = [th2 x(1)] ;
    th3 = [th3 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end

elseif Lmin == r3
    
th2=[];
th1=[];
theta3max = 0;
tmd = rad2deg(theta3max)
tmdi = -360;
time = -(tmd - tmdi)/w
t=0:0.1:time



th3=( w*t + tmd)



%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3);


%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th3)

    theta3 = th3(i)
    F = @(x) [+r1*cosd(x(2)) + r2*cosd(x(1)) - r3*cosd(theta3) - r0;
             +r1*sind(x(2)) + r2*sind(x(1)) - r3*sind(theta3)];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0);
    th2 = [th2 x(1)] ;
    th1 = [th1 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end
   
    
elseif Lmin == r2
    

th2=[];
th3=[];
h = app.enterinitialangularvelocityEditField.Value;
w = rad2deg(-(h));

theta1max = acos((r0^2 + r1^2 -  ((r3 + r2)^2))/(2*r1*r0))

tmd = rad2deg(theta1max)
tmdi = -tmd
time = -(tmd - tmdi)/w
t=0:0.1:time



th1=( w*t + tmd)



%Vector loop equations

F(1) = r1*cos(theta1) + r2*cos(theta2) == r0+ r3*cos(theta3)
F(2) = r1*sin(theta1) + r2*sin(theta2) == 0+ r3*sin(theta3)

%Create the initial point x0.  

x0 = [360 - (180 - (tmd/2));180 - (tmd/2)]

    
for i=1:length(th1)

    theta1 = th1(i)
    F = @(x) [+r1*cosd(theta1) + r2*cosd(x(1)) - r3*cosd(x(2)) - r0;
             +r1*sind(theta1) + r2*sind(x(1)) - r3*sind(x(2))];

     
    %Set options to return iterative display.

    options = optimoptions('fsolve','Display','iter');
    
    %Solve the equations.

    [x,fval] = fsolve(F,x0,options);
    th2 = [th2 x(1)] ;
    th3 = [th3 x(2)] ;
    x0 = [double(x(1));double(x(2))]
    disp( "pop")
end
th1
th2
th3
    
    
end

%Velocity
for i = 1:length(th1)
    A = [ r2*cosd(th2(i)) ,- r3*cosd(th3(i));...
          -r2*sind(th2(i)) ,+ r3*sind(th3(i)) ];
    B = [-r1*w*cosd(th1(i));
        r1*w*sind(th1(i))];
    x = inv(A)*B;
    w2(i) = x(1,1);
    w3(i) = x(2,1);
end

%Acceleration
for i = 1:length(th1)
    A = [ r2*cosd(th2(i)) ,- r3*cosd(th3(i));...
          -r2*sind(th2(i)) ,+ r3*sind(th3(i)) ];
    
    B = [r1*w^2*sind(th1(i)) + r2*w2(i)^2*sind(th2(i)) - r3*w3(i)^2*sind(th3(i));...
        r1*w^2*cosd(th1(i)) + r2*w2(i)^2*cosd(th2(i)) - r3*w3(i)^2*cosd(th3(i))];
    x = A\B;
    a2(i) = x(1,1);
    a3(i) = x(2,1);
end

t1 = t;
t1(1) = [];
w2(1) = [];
w3(1) = [];
a2(1) = [];
a3(1) = [];


function [Lmax,Lmin,L,indx_max,indx_min,indx] = sort1(L0,L1,L2,L3)
A = [L0,L1,L2,L3];
[Lmax,indx_max] = max(A);
[Lmin,indx_min] = min(A);
i=1;
k=1;
while(i<=length(A))
    if(i~=indx_max && i~=indx_min)
        L(k) = A(i);
        indx(k) = i;
        k = k+1;
    end
    i=i+1;
end
end

            plot(app.UIAxes4,t1,a2,'LineWidth',2);
            cla(app.UIAxes5);
            cla(app.UIAxes6);
            cla(app.UIAxes3);
            cla(app.UIAxes2);
            plot(app.UIAxes4_2,t1,a3,'LineWidth',2);
        end

        % Button pushed function: CloseButton
        function CloseButtonPushed(app, event)
            app.UIFigure.Visible = 'off';
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 638 811];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.Scrollable = 'on';

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, 'omega 2 vs time')
            xlabel(app.UIAxes2, 'Time (s)')
            ylabel(app.UIAxes2, 'omega 2 (w2)')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Visible = 'off';
            app.UIAxes2.Position = [14 69 306 185];

            % Create UIAxes3
            app.UIAxes3 = uiaxes(app.UIFigure);
            title(app.UIAxes3, 'omega 3 vs time')
            xlabel(app.UIAxes3, 'time (s)')
            ylabel(app.UIAxes3, 'omega 3 (w3)')
            zlabel(app.UIAxes3, 'Z')
            app.UIAxes3.Visible = 'off';
            app.UIAxes3.Position = [328 69 300 185];

            % Create UIAxes4
            app.UIAxes4 = uiaxes(app.UIFigure);
            title(app.UIAxes4, 'alpha 2 vs time')
            xlabel(app.UIAxes4, 'time (s)')
            ylabel(app.UIAxes4, 'alpha 2 (a2)')
            zlabel(app.UIAxes4, 'Z')
            app.UIAxes4.Visible = 'off';
            app.UIAxes4.Position = [20 69 300 185];

            % Create UIAxes4_2
            app.UIAxes4_2 = uiaxes(app.UIFigure);
            title(app.UIAxes4_2, 'alpha 3 vs time')
            xlabel(app.UIAxes4_2, 'time (s)')
            ylabel(app.UIAxes4_2, 'alpha 3 (a3)')
            zlabel(app.UIAxes4_2, 'Z')
            app.UIAxes4_2.Visible = 'off';
            app.UIAxes4_2.Position = [328 69 300 185];

            % Create UIAxes5
            app.UIAxes5 = uiaxes(app.UIFigure);
            title(app.UIAxes5, 'variation of theta 2 with time')
            xlabel(app.UIAxes5, 'time (s)')
            ylabel(app.UIAxes5, 'theta 2')
            zlabel(app.UIAxes5, 'Z')
            app.UIAxes5.Visible = 'off';
            app.UIAxes5.Position = [20 69 300 185];

            % Create UIAxes6
            app.UIAxes6 = uiaxes(app.UIFigure);
            title(app.UIAxes6, 'variation of theta 3 with time')
            xlabel(app.UIAxes6, 'time (s)')
            ylabel(app.UIAxes6, 'theta 3 ')
            zlabel(app.UIAxes6, 'Z')
            app.UIAxes6.Visible = 'off';
            app.UIAxes6.Position = [328 69 300 185];

            % Create enterlengthofr0EditFieldLabel
            app.enterlengthofr0EditFieldLabel = uilabel(app.UIFigure);
            app.enterlengthofr0EditFieldLabel.HorizontalAlignment = 'center';
            app.enterlengthofr0EditFieldLabel.Position = [48 654 120 22];
            app.enterlengthofr0EditFieldLabel.Text = 'enter length of r0';

            % Create enterlengthofr0EditField
            app.enterlengthofr0EditField = uieditfield(app.UIFigure, 'numeric');
            app.enterlengthofr0EditField.Position = [183 654 100 22];

            % Create enterlengthofr2EditField_2Label
            app.enterlengthofr2EditField_2Label = uilabel(app.UIFigure);
            app.enterlengthofr2EditField_2Label.HorizontalAlignment = 'center';
            app.enterlengthofr2EditField_2Label.Position = [48 583 120 22];
            app.enterlengthofr2EditField_2Label.Text = 'enter length of r2';

            % Create enterlengthofr2EditField_2
            app.enterlengthofr2EditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.enterlengthofr2EditField_2.Position = [183 583 100 22];

            % Create enterlengthofr1EditField_2Label
            app.enterlengthofr1EditField_2Label = uilabel(app.UIFigure);
            app.enterlengthofr1EditField_2Label.HorizontalAlignment = 'center';
            app.enterlengthofr1EditField_2Label.Position = [48 614 120 22];
            app.enterlengthofr1EditField_2Label.Text = 'enter length of r1';

            % Create enterlengthofr1EditField_2
            app.enterlengthofr1EditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.enterlengthofr1EditField_2.Position = [183 614 100 22];

            % Create enterlengthofr3EditField_2Label
            app.enterlengthofr3EditField_2Label = uilabel(app.UIFigure);
            app.enterlengthofr3EditField_2Label.HorizontalAlignment = 'center';
            app.enterlengthofr3EditField_2Label.Position = [48 548 120 22];
            app.enterlengthofr3EditField_2Label.Text = 'enter length of r3';

            % Create enterlengthofr3EditField_2
            app.enterlengthofr3EditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.enterlengthofr3EditField_2.Position = [183 548 100 22];

            % Create checkforchangepointButton
            app.checkforchangepointButton = uibutton(app.UIFigure, 'push');
            app.checkforchangepointButton.ButtonPushedFcn = createCallbackFcn(app, @checkforchangepointButtonPushed, true);
            app.checkforchangepointButton.Position = [32 468 136 22];
            app.checkforchangepointButton.Text = 'check for change point';

            % Create clearButton
            app.clearButton = uibutton(app.UIFigure, 'push');
            app.clearButton.ButtonPushedFcn = createCallbackFcn(app, @clearButtonPushed, true);
            app.clearButton.Position = [208 468 100 22];
            app.clearButton.Text = 'clear';

            % Create SIMULATEButton
            app.SIMULATEButton = uibutton(app.UIFigure, 'push');
            app.SIMULATEButton.ButtonPushedFcn = createCallbackFcn(app, @SIMULATEButtonPushed, true);
            app.SIMULATEButton.FontSize = 16;
            app.SIMULATEButton.FontWeight = 'bold';
            app.SIMULATEButton.Position = [271 410 147 36];
            app.SIMULATEButton.Text = 'SIMULATE ';

            % Create TextArea
            app.TextArea = uitextarea(app.UIFigure);
            app.TextArea.HorizontalAlignment = 'center';
            app.TextArea.Visible = 'off';
            app.TextArea.Position = [415 616 150 60];
            app.TextArea.Value = {'Given mechanism is not change point'};

            % Create TextArea_2
            app.TextArea_2 = uitextarea(app.UIFigure);
            app.TextArea_2.HorizontalAlignment = 'center';
            app.TextArea_2.Visible = 'off';
            app.TextArea_2.Position = [417 616 150 60];
            app.TextArea_2.Value = {'Given mechanism is a change point mechanism '};

            % Create checktypeButton
            app.checktypeButton = uibutton(app.UIFigure, 'push');
            app.checktypeButton.ButtonPushedFcn = createCallbackFcn(app, @checktypeButtonPushed, true);
            app.checktypeButton.Position = [440 574 100 22];
            app.checktypeButton.Text = 'check type';

            % Create TextArea_3
            app.TextArea_3 = uitextarea(app.UIFigure);
            app.TextArea_3.HorizontalAlignment = 'center';
            app.TextArea_3.Visible = 'off';
            app.TextArea_3.Position = [415 489 150 60];
            app.TextArea_3.Value = {'normal model '; '(all side unequal)'};

            % Create TextArea_4
            app.TextArea_4 = uitextarea(app.UIFigure);
            app.TextArea_4.HorizontalAlignment = 'center';
            app.TextArea_4.Visible = 'off';
            app.TextArea_4.Position = [415 489 150 60];
            app.TextArea_4.Value = {'deltoid model'};

            % Create TextArea_5
            app.TextArea_5 = uitextarea(app.UIFigure);
            app.TextArea_5.HorizontalAlignment = 'center';
            app.TextArea_5.Visible = 'off';
            app.TextArea_5.Position = [417 489 150 60];
            app.TextArea_5.Value = {'parallelogram model'};

            % Create ANALYSISButton
            app.ANALYSISButton = uibutton(app.UIFigure, 'push');
            app.ANALYSISButton.ButtonPushedFcn = createCallbackFcn(app, @ANALYSISButtonPushed, true);
            app.ANALYSISButton.FontSize = 16;
            app.ANALYSISButton.FontWeight = 'bold';
            app.ANALYSISButton.Position = [274 343 141 40];
            app.ANALYSISButton.Text = 'ANALYSIS ';

            % Create TextArea_6
            app.TextArea_6 = uitextarea(app.UIFigure);
            app.TextArea_6.HorizontalAlignment = 'center';
            app.TextArea_6.Visible = 'off';
            app.TextArea_6.Position = [417 489 150 60];
            app.TextArea_6.Value = {'parallelogram '; '(all side equal)'};

            % Create velocityanalysisButton
            app.velocityanalysisButton = uibutton(app.UIFigure, 'push');
            app.velocityanalysisButton.ButtonPushedFcn = createCallbackFcn(app, @velocityanalysisButtonPushed, true);
            app.velocityanalysisButton.Position = [270 291 131 22];
            app.velocityanalysisButton.Text = 'velocity analysis ';

            % Create accelerationanalysisButton_2
            app.accelerationanalysisButton_2 = uibutton(app.UIFigure, 'push');
            app.accelerationanalysisButton_2.ButtonPushedFcn = createCallbackFcn(app, @accelerationanalysisButton_2Pushed, true);
            app.accelerationanalysisButton_2.Position = [434 291 131 22];
            app.accelerationanalysisButton_2.Text = 'acceleration analysis ';

            % Create CloseButton
            app.CloseButton = uibutton(app.UIFigure, 'push');
            app.CloseButton.ButtonPushedFcn = createCallbackFcn(app, @CloseButtonPushed, true);
            app.CloseButton.Position = [270 26 100 22];
            app.CloseButton.Text = 'Close';

            % Create GROUP14Label
            app.GROUP14Label = uilabel(app.UIFigure);
            app.GROUP14Label.FontSize = 16;
            app.GROUP14Label.FontWeight = 'bold';
            app.GROUP14Label.Position = [276 706 87 22];
            app.GROUP14Label.Text = 'GROUP-14';

            % Create GUIFORCHANGEPOINTMECHANISMANALYSISLabel
            app.GUIFORCHANGEPOINTMECHANISMANALYSISLabel = uilabel(app.UIFigure);
            app.GUIFORCHANGEPOINTMECHANISMANALYSISLabel.FontSize = 16;
            app.GUIFORCHANGEPOINTMECHANISMANALYSISLabel.FontWeight = 'bold';
            app.GUIFORCHANGEPOINTMECHANISMANALYSISLabel.Position = [127 739 397 22];
            app.GUIFORCHANGEPOINTMECHANISMANALYSISLabel.Text = 'GUI FOR CHANGE POINT MECHANISM ANALYSIS ';

            % Create PHY113Label
            app.PHY113Label = uilabel(app.UIFigure);
            app.PHY113Label.FontSize = 16;
            app.PHY113Label.FontWeight = 'bold';
            app.PHY113Label.Position = [276 727 387 113];
            app.PHY113Label.Text = '19PHY113';

            % Create enterinitialangularvelocityEditFieldLabel
            app.enterinitialangularvelocityEditFieldLabel = uilabel(app.UIFigure);
            app.enterinitialangularvelocityEditFieldLabel.HorizontalAlignment = 'right';
            app.enterinitialangularvelocityEditFieldLabel.Position = [14 508 154 22];
            app.enterinitialangularvelocityEditFieldLabel.Text = 'enter initial angular velocity ';

            % Create enterinitialangularvelocityEditField
            app.enterinitialangularvelocityEditField = uieditfield(app.UIFigure, 'numeric');
            app.enterinitialangularvelocityEditField.Position = [183 508 100 22];

            % Create thetaButton
            app.thetaButton = uibutton(app.UIFigure, 'push');
            app.thetaButton.ButtonPushedFcn = createCallbackFcn(app, @thetaButtonPushed, true);
            app.thetaButton.Position = [109 291 100 22];
            app.thetaButton.Text = 'theta';

            % Create Label
            app.Label = uilabel(app.UIFigure);
            app.Label.Position = [10 685 651 22];
            app.Label.Text = '--------------------------------------------------------------------------------------------------------------------------------------------------------------';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = mech_analyser_exported_2

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end