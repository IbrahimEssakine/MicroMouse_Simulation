function varargout = InterfaceMazeSolver(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @InterfaceMazeSolver_OpeningFcn, ...
                   'gui_OutputFcn',  @InterfaceMazeSolver_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
end
% End initialization code - DO NOT EDIT


% --- Executes just before InterfaceMazeSolver is made visible.
function InterfaceMazeSolver_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InterfaceMazeSolver (see VARARGIN)

% Choose default command line output for InterfaceMazeSolver
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
axes(handles.FullMaze);
set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
axes(handles.MapRobot);
set(gca, 'XTick', [], 'YTick', []); % Hide axis labels

I = imread('/2.jpg');
axes(handles.FullMaze);
imshow(I);
I = imread('/1.jpg');
axes(handles.MapRobot);
imshow(I);
 readymazematrix =[];
handles.readymazematrix=readymazematrix;
guidata(hObject,handles)

sz =0;
handles.sz=sz;
guidata(hObject,handles)

MazeComplet =[];
handles.MazeComplet=MazeComplet;
guidata(hObject,handles)

end

% UIWAIT makes InterfaceMazeSolver wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = InterfaceMazeSolver_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end


% --- Executes on button press in Solve.
function Solve_Callback(hObject, eventdata, handles)
% hObject    handle to Solve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  readymazematrix=handles.readymazematrix;
    maze=readymazematrix ;
%  axes(handles.MapRobot);
 
 
    sz=handles.sz;
    mid=fix(sz/2);
    directions = [0, 1; 1, 0; 0, -1; -1, 0];
%     choose starting point 
    [sy,sx]=ginput(1);
    sx=sz-fix(sx)+1;
    sy=fix(sy);
    maze(sx,sy)=4;
    
    %addition to mouad's idea
    MazeBlind=ones(sz,sz);
    MazeBlind(sx,sy)=4;
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.MapRobot,FlipedMazeBlind);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    
                                    
%     choose ending point
    [ey,ex]=ginput(1);
    ex=sz-fix(ex)+1;
    ey=fix(ey);
%     2 is the identifier of the ending point = destination
    maze(ex,ey)=2;
    
    handles.MazeComplet=maze;
    MazeBlind(sx,sy)=2;

    queue=[sx,sy];
    exception=0;
    parent = cell(sz,sz);
    tic
    while size(queue)>0 

        for k=1:size(queue,1)
            x=queue(1,1);
            y=queue(1,2);
            queue=queue(2:end,:);
            neighbors = [];
            
            for i = 1:4
                nx = x + directions(i, 1);
                ny = y + directions(i, 2);
                if nx > 0 && ny > 0 && nx <= sz && ny <= sz && maze(nx,ny) ~=3
                    neighbors = [neighbors; nx, ny,i];
                end
            end

             if  ~isempty(neighbors)
                for i = 1:size(neighbors)
                    neighbor = neighbors(i, :);
                    nx = neighbor(1);
                    ny = neighbor(2);
                    
                    MazeBlind(nx,ny)=maze(nx,ny);
                    if (maze(nx, ny) == 1)
                        MazeBlind(nx,ny)=3;
                        maze(nx,ny)=3;
                    end

                    if maze(nx, ny) == 0
                        continue
                    end
                    maze(nx,ny)
                    if maze(nx,ny)==2
                        destination = neighbor;
                        parent{nx, ny} = [x,y];
                        exception =1;
                        break
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    FlipedMazeBlind = flipud(MazeBlind);
                    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
                    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
                    animation=findobj('Tag', 'animation');
                    if (animation.Value==1)
                        pcolor(handles.MapRobot,FlipedMazeBlind);
                    end
                    
                    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
                    pause(0);
                    %%%%%%%%%%%%%%%%%%
                    
                    parent{nx, ny} = [x,y];
                    
                    queue= [queue; nx,ny];
                end
             end
             if exception ==1
                 break
             end


        end
        if exception ==1
            break
        end
    end
    timer=toc;
    TimerText = findobj('Tag', 'Exploration');
    set(TimerText, 'String', num2str(timer));
    
    
    path=destination;
    pathValue=1;
    tic
    while (path(end,1)~=sx || path(end,2)~=sy) 
        pathValue=pathValue+1;
        cmap = [

            0.2,0.2,0.2;
            1, 1, 1;
            1,1,0;
            0,1,0;
            0.1,0.7,0.5;
        ];
        colormap(cmap);

        FlipedMazeBlind = flipud(MazeBlind);
        FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
        FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
        animation=findobj('Tag', 'animation');
        if (animation.Value==1)
            pcolor(handles.MapRobot,FlipedMazeBlind);
        end
        
        set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
        pause(0)
        
        
       cellparent=parent{path(end,1),path(end,2)};
       cx=cellparent(end,1);
       cy=cellparent(end,2);
       MazeBlind(cx,cy)=4;
         maze(cx,cy)=4;
       path=cellparent; 
    end
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.MapRobot,FlipedMazeBlind);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    timer=toc;
    ScoreText = findobj('Tag', 'Score');
    set(ScoreText, 'String', num2str(pathValue));
    TimerText = findobj('Tag', 'Solution');
    set(TimerText, 'String', num2str(timer));
      FlipedMazeBlind = flipud(maze);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.FullMaze,FlipedMazeBlind);
end
% --- Executes on button press in LoadMap.
function LoadMap_Callback(hObject, eventdata, handles)
% hObject    handle to LoadMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 [filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp;*.tiff;*.tif;*.gif', 'Image Files (*.jpg, *.png, *.bmp, *.tiff, *.tif, *.gif)'; '*.*', 'All Files (*.*)'}, 'Select Image');

    if isequal(filename,0) || isequal(pathname,0)
       disp('User pressed cancel')
    else
       imagefile = strcat(pathname,filename);
       I = imread(imagefile);
       axes(handles.FullMaze);
       imshow(I);
       handles.I = I;
       guidata(hObject, handles);
    end
    
    sz=str2double(inputdlg('Entre The Resolution : '));
    handles.sz=sz;
    img= imread(imagefile);

    img = rgb2gray(img);

    % Binarize the image
    bin_img = imbinarize(img);

    % If you want 0s to represent white and 1s to represent black, you can invert the binary image
    % bin_img = ~bin_img;


    % This is done by taking the mean of each block of pixels that represents a cell of the maze
    % The round function is used to convert the mean values to 0 or 1
    mat = round(imresize(bin_img, [sz sz], 'bicubic'));
    
%     mid=fix(sz/2);
%     mat(mid+1+[0,1],mid+1+[0,1])=2;
    cmap = [

    0.2,0.2,0.2;
    1, 1, 1;
    1,1,0;
    0,1,0;
];
     colormap(cmap);
    
    % Flip the maze matrix in the up-down direction
    maze = flipud(mat);
    maze(end+1,:) = maze(end,:); % Append last row
    maze(:,end+1) = maze(:,end); % Append last column
    % Create a `pcolor` plot of the matrix

     readymazematrix = mat;
    handles.readymazematrix=readymazematrix;
    guidata(hObject,handles)
    
    pcolor(handles.FullMaze , maze);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels

end
   
 
% --- Executes on button press in SolveDFS.
function SolveDFS_Callback(hObject, eventdata, handles)
readymazematrix=handles.readymazematrix;
    maze=readymazematrix ;
%  axes(handles.MapRobot);
 
 
    sz=handles.sz;
    mid=fix(sz/2);
    directions = [0, 1; 1, 0; 0, -1; -1, 0];
%     choose starting point 
    [sy,sx]=ginput(1);
    sx=sz-fix(sx)+1;
    sy=fix(sy);
    maze(sx,sy)=4;
    
    %addition to mouad's idea
    MazeBlind=ones(sz,sz);
    MazeBlind(sx,sy)=4;
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.MapRobot,FlipedMazeBlind);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    
                                    
%     choose ending point

    [ey,ex]=ginput(1);
    ex=sz-fix(ex)+1;
    ey=fix(ey);
%   2 is the identifier of the ending point = destination
    maze(ex,ey)=2;
    handles.readymazematrix=maze;
    MazeBlind(sx,sy)=2;
    instantDir=[1,0];
    previousPoint=[sx1,sy];
    queue=[sx,sy,sz*sz*virageSpeed,speed,instantDir,previousPoint];
    exception=0;
    parent = cell(sz,sz);
    tic
    while size(queue)>0 

        for k=1:size(queue,1)
            x=queue(1,1);
            y=queue(1,2);
            timePassed=queue(1,3);
            instantSpeed=queue(1,4);
            instantDir=queue(1,5);
            queue=queue(2:end,:);
            neighbors = [];
            for i = 1:4
                nx = x + directions(i, 1);
                ny = y + directions(i, 2);
                if nx > 0 && ny > 0 && nx <= sz && ny <= sz && maze(nx,ny) ~=3
                    neighbors = [neighbors; nx, ny,i];
                end
            end
            
             if  ~isempty(neighbors)
                for i = 1:size(neighbors)
                    neighbor = neighbors(i, :);
                    nx = neighbor(1);
                    ny = neighbor(2);
                    new_direction = updateDirection(previous_point, new_point, previous_direction);
                    MazeBlind(nx,ny)=maze(nx,ny);
                    if (maze(nx, ny) == 1)
                        MazeBlind(nx,ny)=3;
                        maze(nx,ny)=3;
                    end

                    if maze(nx, ny) == 0
                        continue
                    end
                    if maze(nx,ny)==2
                        destination = neighbor;
                        parent{nx, ny} = [x,y];
                        exception =1;
                        break
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    FlipedMazeBlind = flipud(MazeBlind);
                    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
                    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
                    animation=findobj('Tag', 'animation');
                    if (animation.Value==1)
                        pcolor(handles.MapRobot,FlipedMazeBlind);
                    end
                    
                    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
                    pause(0);
                    %%%%%%%%%%%%%%%%%%
                    
                    parent{nx, ny} = [x,y];
                    
                    queue=[nx,ny,timePassed,instantSpeed,instantDir];

                    % Sort the array based on the 3rd column (3rd element) in ascending order
                    queue = sortrows(queue, -3);
                end
             end
             if exception ==1
                 break
             end


        end
        if exception ==1
            break
        end
    end
    pcolor(handles.MapRobot,FlipedMazeBlind);
    timer=toc;
    TimerText = findobj('Tag', 'Exploration');
    set(TimerText, 'String', num2str(timer));
    
    
    path=destination;
    pathValue=1;
    tic
    while (path(end,1)~=sx || path(end,2)~=sy) 
        pathValue=pathValue+1;
        cmap = [

            0.2,0.2,0.2;
            1, 1, 1;
            1,1,0;
            0,1,0;
            1,0.1,0.1;
        ];
        colormap(cmap);

        FlipedMazeBlind = flipud(MazeBlind);
        FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
        FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
        animation=findobj('Tag', 'animation');
        if (animation.Value==1)
            pcolor(handles.MapRobot,FlipedMazeBlind);
        end
        set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
        pause(0);
        
        
       cellparent=parent{path(end,1),path(end,2)};
       cx=cellparent(end,1);
       cy=cellparent(end,2);
       MazeBlind(cx,cy)=4;
         maze(cx,cy)=4;
       path=cellparent; 
    end
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.MapRobot,FlipedMazeBlind);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    timer=toc;
    ScoreText = findobj('Tag', 'Score');
    set(ScoreText, 'String', num2str(pathValue));
    TimerText = findobj('Tag', 'Solution');
    set(TimerText, 'String', num2str(timer));
    FlipedMazeBlind = flipud(maze);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.FullMaze,FlipedMazeBlind);

end

% --- Executes on button press in BreakWalls.
function BreakWalls_Callback(hObject, eventdata, handles)
% hObject    handle to BreakWalls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    % hObject    handle to Solve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  readymazematrix=handles.readymazematrix;
    maze=readymazematrix ;
%  axes(handles.MapRobot);
 
    secMaze=maze;
    sz=handles.sz;
    mid=fix(sz/2);
    directions = [0, 1; 1, 0; 0, -1; -1, 0];
%     choose starting point 
    [sy,sx]=ginput(1);
    sx=sz-fix(sx)+1;
    sy=fix(sy);
    maze(sx,sy)=4;
    
    %addition to mouad's idea
    MazeBlind=ones(sz,sz);
    MazeBlind(sx,sy)=4;
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.MapRobot,FlipedMazeBlind);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    secMaze=MazeBlind;
                                    
%     choose ending point

    [ey,ex]=ginput(1);
    ex=sz-fix(ex)+1;
    ey=fix(ey);
%     2 is the identifier of the ending point = destination
    maze(ex,ey)=2;
        
    handles.MazeComplet=maze;
    MazeBlind(sx,sy)=2;


   
        rows = 2:size(maze,1) - 1;
        cols = 2:size(maze,1) - 1;
        sz=size(maze,1);
        % Initialize a flag to track whether a zero was changed

        % Loop until a zero is changed
        
        zeroIndices = find(maze(rows, cols) == 0);
        values=[];
        for i=1:size((zeroIndices),1)
            incideMaze=maze(rows,cols);
            incideMaze(zeroIndices(i))=1;
            secMaze=maze;
            secMaze(rows,cols)=incideMaze;
            zeroIndice=zeroIndices(i);
            
            
            queue=[sx,sy];
    exception=0;
    parent = cell(sz,sz);
    tic
    while size(queue)>0 

        for k=1:size(queue,1)
            x=queue(1,1);
            y=queue(1,2);
            queue=queue(2:end,:);
            neighbors = [];
            
            for i = 1:4
                nx = x + directions(i, 1);
                ny = y + directions(i, 2);
                if nx > 0 && ny > 0 && nx <= sz && ny <= sz && secMaze(nx,ny) ~=3
                    neighbors = [neighbors; nx, ny,i];
                end
            end

             if  ~isempty(neighbors)
                for i = 1:size(neighbors)
                    neighbor = neighbors(i, :);
                    nx = neighbor(1);
                    ny = neighbor(2);
                    
                    MazeBlind(nx,ny)=secMaze(nx,ny);
                    if (secMaze(nx, ny) == 1)
                        MazeBlind(nx,ny)=3;
                        secMaze(nx,ny)=3;
                    end

                    if secMaze(nx, ny) == 0
                        continue
                    end
                    if secMaze(nx,ny)==2
                        destination = neighbor;
                        parent{nx, ny} = [x,y];
                        exception =1;
                        break
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    FlipedMazeBlind = flipud(MazeBlind);
                    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
                    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
                    animation=findobj('Tag', 'animation');
                    if (animation.Value==1)
                        pcolor(handles.MapRobot,FlipedMazeBlind);
                    end
                    
                    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
                    pause(0);
                    %%%%%%%%%%%%%%%%%%
                    
                    parent{nx, ny} = [x,y];
                    
                    queue= [queue; nx,ny];
                end
             end
             if exception ==1
                 break
             end


        end
        if exception ==1
            break
        end
    end
    timer=toc;
    TimerText = findobj('Tag', 'Exploration');
    set(TimerText, 'String', num2str(timer));
    
    
    path=destination;
    pathValue=1;
    tic
    while (path(end,1)~=sx || path(end,2)~=sy) 
        pathValue=pathValue+1;
        cmap = [

            0.2,0.2,0.2;
            1, 1, 1;
            1,1,0;
            0,1,0;
            0.1,0.7,0.5;
        ];
        colormap(cmap);

        FlipedMazeBlind = flipud(MazeBlind);
        FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
        FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
        animation=findobj('Tag', 'animation');
        if (animation.Value==1)
            pcolor(handles.MapRobot,FlipedMazeBlind);
        end
        
        set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
        pause(0)
        
        
       cellparent=parent{path(end,1),path(end,2)};
       cx=cellparent(end,1);
       cy=cellparent(end,2);
       MazeBlind(cx,cy)=4;
       secMaze(cx,cy)=4;
       path=cellparent; 
    end
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.MapRobot,FlipedMazeBlind);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    timer=toc;
    ScoreText = findobj('Tag', 'Score');
    set(ScoreText, 'String', num2str(pathValue));
    TimerText = findobj('Tag', 'Solution');
    set(TimerText, 'String', num2str(timer));
      FlipedMazeBlind = flipud(secMaze);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.FullMaze,FlipedMazeBlind);
            
            values = [values;pathValue zeroIndice];
            if (abs(sx-ex)+abs(sy-ey))==(pathValue-1)
                break;
            end
        end
        SortedValues=sortrows(values, 1);
        breakPoint=SortedValues(1,2);
        
        incideMaze=maze(rows,cols);
        incideMaze(breakPoint)=1;
        secMaze=maze;
        secMaze(rows,cols)=incideMaze;
        queue=[sx,sy];
    exception=0;
    parent = cell(sz,sz);
    tic
    while size(queue)>0 

        for k=1:size(queue,1)
            x=queue(1,1);
            y=queue(1,2);
            queue=queue(2:end,:);
            neighbors = [];
            
            for i = 1:4
                nx = x + directions(i, 1);
                ny = y + directions(i, 2);
                if nx > 0 && ny > 0 && nx <= sz && ny <= sz && secMaze(nx,ny) ~=3
                    neighbors = [neighbors; nx, ny,i];
                end
            end

             if  ~isempty(neighbors)
                for i = 1:size(neighbors)
                    neighbor = neighbors(i, :);
                    nx = neighbor(1);
                    ny = neighbor(2);
                    
                    MazeBlind(nx,ny)=secMaze(nx,ny);
                    if (secMaze(nx, ny) == 1)
                        MazeBlind(nx,ny)=3;
                        secMaze(nx,ny)=3;
                    end

                    if secMaze(nx, ny) == 0
                        continue
                    end
                    if secMaze(nx,ny)==2
                        destination = neighbor;
                        parent{nx, ny} = [x,y];
                        exception =1;
                        break
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    FlipedMazeBlind = flipud(MazeBlind);
                    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
                    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
                    animation=findobj('Tag', 'animation');
                    if (animation.Value==1)
                        pcolor(handles.MapRobot,FlipedMazeBlind);
                    end
                    
                    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
                    pause(0);
                    %%%%%%%%%%%%%%%%%%
                    
                    parent{nx, ny} = [x,y];
                    
                    queue= [queue; nx,ny];
                end
             end
             if exception ==1
                 break
             end


        end
        if exception ==1
            break
        end
    end
    timer=toc;
    TimerText = findobj('Tag', 'Exploration');
    set(TimerText, 'String', num2str(timer));
    
    
    path=destination;
    pathValue=1;
    tic
    while (path(end,1)~=sx || path(end,2)~=sy) 
        pathValue=pathValue+1;
        cmap = [

            0.2,0.2,0.2;
            1, 1, 1;
            1,1,0;
            0,1,0;
            0.1,0.7,0.5;
        ];
        colormap(cmap);

        FlipedMazeBlind = flipud(MazeBlind);
        FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
        FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
        animation=findobj('Tag', 'animation');
        if (animation.Value==1)
            pcolor(handles.MapRobot,FlipedMazeBlind);
        end
        
        set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
        pause(0)
        
        
       cellparent=parent{path(end,1),path(end,2)};
       cx=cellparent(end,1);
       cy=cellparent(end,2);
       MazeBlind(cx,cy)=4;
       secMaze(cx,cy)=4;
       path=cellparent; 
    end
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.MapRobot,FlipedMazeBlind);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    timer=toc;
    ScoreText = findobj('Tag', 'Score');
    set(ScoreText, 'String', num2str(pathValue));
    TimerText = findobj('Tag', 'Solution');
    set(TimerText, 'String', num2str(timer));
      FlipedMazeBlind = flipud(secMaze);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.FullMaze,FlipedMazeBlind);
    
            
end

% --- Executes on button press in AStar.

function AStar_Callback(hObject, eventdata, handles)
% hObject    handle to AStar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

readymazematrix=handles.readymazematrix;
    maze=readymazematrix ;
%  axes(handles.MapRobot);
 
 
    sz=handles.sz;
    mid=fix(sz/2);
    directions = [0, 1; 1, 0; 0, -1; -1, 0];
%     choose starting point 
    [sy,sx]=ginput(1);
    sx=sz-fix(sx)+1;
    sy=fix(sy);
    maze(sx,sy)=4;
    
    %addition to mouad's idea
    MazeBlind=ones(sz,sz);
    MazeBlind(sx,sy)=4;
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.MapRobot,FlipedMazeBlind);
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    
    
   
                                    
%     choose ending point

    [ey,ex]=ginput(1);
    ex=sz-fix(ex)+1;
    ey=fix(ey);
%     2 is the identifier of the ending point = destination
    maze(ex,ey)=2;
    handles.readymazematrix=maze;
    MazeBlind(sx,sy)=2;
    % making the g_score and f_score
    g_score=inf(sz);
    g_score(sx,sy)=0;
    h_score=inf(sz);
    for i=1:sz
        for j=1:sz
            if (maze(i,j)== 0)
                continue;
            end
            h_score(i,j)=abs(i-ex)+abs(j-ey);
        end
    end
    f_score=inf(sz);
    f_score(sx,sy)=h_score(sx,sy);
    queue={abs(sx-0)+abs(sy-0)+f_score(sx,sy) h_score(sx,sy) sx sy;inf inf 0 0};
    %currentPosition={sx,sy};
    exception=0;
    parent = cell(sz,sz);
    last_x=0;
    last_y=0;
    tic
    while size(queue)>0 
            %current Position
            x=cell2mat(queue(1,3));
            y=cell2mat(queue(1,4));
            queue(1,:)=[];
            if(last_x==x && last_y==y)
                queue(1, 1) = {inf};
                queue = sortrows(queue, [1, 2]);
                continue;
            end
            last_x=x;
            last_y=y;
            
            neighbors = [];
            for i = 1:4
                nx = x + directions(i, 1);
                ny = y + directions(i, 2);
                MazeBlind(nx,ny)=maze(nx,ny);
                if nx > 0 && ny > 0 && nx <= sz && ny <= sz && maze(nx,ny) ~=3 && maze(nx,ny)~=0
                    neighbors = [neighbors; nx, ny,i];
                    temp_g_score=g_score(x,y)+1;
                    temp_f_score=temp_g_score+h_score(nx,ny);
                    if temp_f_score<f_score(nx,ny)
                        g_score(nx,ny)=temp_g_score;
                        f_score(nx,ny)=temp_f_score;
                        queue=[queue;{abs(nx-sx)+abs(sy-ny)+f_score(nx,ny) h_score(nx,ny) nx ny}];
                        queue = sortrows(queue, [1, 2]);
                        
                    end
                        
                end
            end
            
             if  ~isempty(neighbors)
                for i = 1:size(neighbors)
                    neighbor = neighbors(i, :);
                    nx = neighbor(1);
                    ny = neighbor(2);
                    MazeBlind(nx,ny)=maze(nx,ny);
                    if (maze(nx, ny) == 1)
                        MazeBlind(nx,ny)=3;
                        maze(nx,ny)=3;
                    end

                    if maze(nx, ny) == 0
                        continue
                    end
                    
                    if maze(nx,ny)==2
                        destination = neighbor;
                        parent{nx, ny} = [x,y];
                        exception =1;
                        break
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    FlipedMazeBlind = flipud(MazeBlind);
                    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
                    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
                    animation=findobj('Tag', 'animation');
                    if (animation.Value==1)
                        pcolor(handles.MapRobot,FlipedMazeBlind);
                    end
                    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
                    pause(0);
                    %%%%%%%%%%%%%%%%%%
                    
                    parent{nx, ny} = [x,y];
                    

                end
             end
             if exception ==1
                 break
             end


        
        if exception ==1
            break
        end
    end
    timer=toc;
    TimerText = findobj('Tag', 'Exploration');
    set(TimerText, 'String', num2str(timer));
    
    
    path=destination;
    pathValue=1;
    tic
    while (path(end,1)~=sx || path(end,2)~=sy) 
        pathValue=pathValue+1;
        cmap = [

            0.2,0.2,0.2;
            1, 1, 1;
            1,1,0;
            0,1,0;
            1,0.1,0.1;
        ];
        colormap(cmap);

        FlipedMazeBlind = flipud(MazeBlind);
        FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
        FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
        animation=findobj('Tag', 'animation');
        if (animation.Value==1)
            pcolor(handles.MapRobot,FlipedMazeBlind);
        end
        set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
        pause(0);
        
        
       cellparent=parent{path(end,1),path(end,2)};
       cx=cellparent(end,1);
       cy=cellparent(end,2);
       MazeBlind(cx,cy)=4;
       maze(cx,cy)=4;
       path=cellparent; 
    end
    FlipedMazeBlind = flipud(MazeBlind);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    animation=findobj('Tag', 'animation');
    pcolor(handles.MapRobot,FlipedMazeBlind);
    
    
    set(gca, 'XTick', [], 'YTick', []); % Hide axis labels
    timer=toc;
    ScoreText = findobj('Tag', 'Score');
    set(ScoreText, 'String', num2str(pathValue));
    TimerText = findobj('Tag', 'Solution');
    set(TimerText, 'String', num2str(timer));
    FlipedMazeBlind = flipud(maze);
    FlipedMazeBlind(end+1,:) = FlipedMazeBlind(end,:); % Append last row
    FlipedMazeBlind(:,end+1) = FlipedMazeBlind(:,end); % Append last column
    pcolor(handles.FullMaze,FlipedMazeBlind);

end


% --- Executes on button press in Pause_Mode.
function Pause_Mode_Callback(hObject, eventdata, handles)
       pause();
end


% --- Executes on button press in animation.
function animation_Callback(hObject, eventdata, handles)
% hObject    handle to animation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of animation
end