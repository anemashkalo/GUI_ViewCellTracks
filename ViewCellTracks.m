function varargout = ViewCellTracks(varargin)
% VIEWCELLTRACKS MATLAB code for ViewCellTracks.fig
%      VIEWCELLTRACKS, by itself, creates a new VIEWCELLTRACKS or raises the existing
%      singleton*.
%
%      H = VIEWCELLTRACKS returns the handle to a new VIEWCELLTRACKS or the handle to
%      the existing singleton*.
%
%      VIEWCELLTRACKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWCELLTRACKS.M with the given input arguments.
%
%      VIEWCELLTRACKS('Property','Value',...) creates a new VIEWCELLTRACKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewCellTracks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewCellTracks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewCellTracks

% Last Modified by GUIDE v2.5 21-Feb-2018 14:30:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ViewCellTracks_OpeningFcn, ...
                   'gui_OutputFcn',  @ViewCellTracks_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before ViewCellTracks is made visible.
function ViewCellTracks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViewCellTracks (see VARARGIN)

% Choose default command line output for ViewCellTracks
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ViewCellTracks wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ViewCellTracks_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.currT = int32(get(handles.slider1, 'Value'));%

if isempty(handles.tracktoplot) || isempty(handles.tracks)
hold on;UpdateImage(handles,handles.currT);
title(['Current time point  ' num2str((double(handles.currT)*handles.delta_t)/60)  ...
    'hrs since start; Frame(' num2str(handles.currT) ')']);hold on
else
xyt = handles.tracks;    
handles.currT = int32(get(handles.slider1, 'Value'));%
showfulltrackbutton = get(handles.ShowFullTrack,'String');
showtrackIDsbutton = get(handles.ShowTrackIDs,'String');
if showtrackIDsbutton(1) == 'H' % if the button was not depressed before clicking the slider
handles.ShowTrackIDsState = handles.ShowTrackIDsState+1;
set(handles.ShowTrackIDs,'String','ShowTrackIDs');
end
% below part allows the full track to be on, while scrolling, untill the user
% presses hidefulltracks
if showfulltrackbutton(1) == 'H' && ~isempty(handles.tracktoplot)% if the button was not depressed before%clicking the slider%%%%%
toview = size(handles.tracktoplot,2);
xyt_curr = struct;
for jj=1:toview
   xyt_curr(jj).fulltrack = xyt(handles.tracktoplot(jj)).dat;  
end
handles.fulltrack = xyt_curr; 
handles.currT = int32(get(handles.slider1, 'Value'));
UpdateImage(handles,handles.currT);hold on
colormap = handles.colormap;
for k=1:toview
 firstT(k)= min(nonzeros(xyt_curr(k).fulltrack(:,3))); 
 lastT(k) =  max(nonzeros(xyt_curr(k).fulltrack(:,3)));
if (firstT(k)<=handles.currT) && (handles.currT<=lastT(k)) % if this track has an XY at this time point
hplot2 = plot(xyt_curr(k).fulltrack(firstT(k)+1:end-1,1),xyt_curr(k).fulltrack(firstT(k)+1:end-1,2)...
    ,'-*','Color',colormap(handles.trackcolor(k),:),'MarkerSize',3,'LineWidth',1);hold on
plot(xyt_curr(k).fulltrack(firstT(k),1),xyt_curr(k).fulltrack(firstT(k),2)...
    ,'p','MarkerFaceColor','r','MarkerEdgeColor','w','MarkerSize',6,'LineWidth',1);hold on
text(xyt_curr(k).fulltrack(firstT(k),1)+8,xyt_curr(k).fulltrack(firstT(k),2)+5,'Start','Color','r','FontSize',8);hold on
plot(xyt_curr(k).fulltrack(end,1),xyt_curr(k).fulltrack(end,2)...
    ,'p','MarkerFaceColor','g','MarkerEdgeColor','w','MarkerSize',6,'LineWidth',1);hold on
text(xyt_curr(k).fulltrack(end,1)+8,xyt_curr(k).fulltrack(end,2)+5,'End','Color','g','FontSize',8);hold on
% label the current point
if  handles.currT<=size(xyt_curr(k).fulltrack,1)% not all cells are tracked all the way
        if xyt_curr(k).fulltrack(handles.currT,1)>0 % tracks that were assigned after first time point, have zeros in the beginning so don't plot them
            plot(xyt_curr(k).fulltrack(handles.currT,1),xyt_curr(k).fulltrack(handles.currT,2),'p',...
                'MarkerFaceColor',colormap(handles.trackcolor(k),:),'MarkerEdgeColor','k','Markersize',15);hold on  
            text(handles.axes1,xyt_curr(k).fulltrack(handles.currT,1)+5,xyt_curr(k).fulltrack(handles.currT,2)+5,...
                num2str(handles.tracktoplot(k)),'Color',colormap(handles.trackcolor(k),:),'FontSize',12);hold on
        else
            disp(['Track ID ' num2str(handles.tracktoplot(k)) 'was picked up in frame' num2str(firstT(k))]);
        end
end
end
end
title(['Current time point  ' num2str((double(handles.currT)*handles.delta_t)/60)  ...
    'hrs since start; Frame(' num2str(handles.currT) ')']);hold on
else
toview = size(handles.tracktoplot,2);
hold on;UpdateImage(handles,handles.currT);%
xyt_curr = struct;sz_tmp = [];
for jj=1:toview
   xyt_curr(jj).fulltrack = xyt(handles.tracktoplot(jj)).dat;
   sz_tmp(jj) = size(xyt(handles.tracktoplot(jj)).dat,1);   
end
handles.fulltrack = xyt_curr;
colormap = handles.colormap;
for k =1:toview
 firstT(k)= min(nonzeros(xyt_curr(k).fulltrack(:,3)));
    if  handles.currT<=size(xyt_curr(k).fulltrack,1)% not all cells are tracked all the way
        if xyt_curr(k).fulltrack(handles.currT,1)>0 % tracks that were assigned after first time point, have zeros in the beginning so don't plot them
            plot(xyt_curr(k).fulltrack(handles.currT,1),xyt_curr(k).fulltrack(handles.currT,2),'p',...
                'MarkerFaceColor',colormap(handles.trackcolor(k),:),'MarkerEdgeColor',colormap(handles.trackcolor(k),:),'Markersize',5);hold on  
            text(handles.axes1,xyt_curr(k).fulltrack(handles.currT,1)+5,xyt_curr(k).fulltrack(handles.currT,2)+5,...
                num2str(handles.tracktoplot(k)),'Color',colormap(handles.trackcolor(k),:),'FontSize',12);hold on
        else
            disp(['Track ID ' num2str(handles.tracktoplot(k)) 'was picked up in frame' num2str(firstT(k))]);
        end
end
end
title(['Current time point  ' num2str((double(handles.currT)*handles.delta_t)/60)  ...
    'hrs since start; Frame(' num2str(handles.currT) ')']);hold on
end

% do if the checkbox was selected and the matfile is loaded;
if ~isempty(handles.matfile) && (handles.show_qFluorData==1)
    handles.PlotFluorData_SelectedCells.Visible = 'On';
    if isempty(handles.displayedFluorData)
    handles.YlabelString.Visible = 'on';
    set(handles.YlabelString,'String','PlotFluor Y label');
    end
    plotfluor = zeros(size(handles.quantifyselected,2),3);
    colormap2 = hot;% for the florIntensitydata
    for k=1:size(handles.quantifyselected,2)
        if size(handles.quantifyselected(k).selected.dat,2)>=handles.currT % if the track exists at the slider position
        plotfluor(k,1:2) = round(handles.quantifyselected(k).selected.dat(handles.currT).xy);
        plotfluor(k,3) = handles.quantifyselected(k).selected.dat(handles.currT).fluor;
        text(handles.axes1,plotfluor(k,1)+10,plotfluor(k,2)-25,num2str(plotfluor(k,3),2),'Color',colormap(handles.trackcolor(k),:),'FontSize',12);hold on  
        end              
    end    
end
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
handles.ShowFullTrack.Visible = 'On';


% --- Executes on button press in loadmaxprojections.
function loadmaxprojections_Callback(hObject, eventdata, handles)
% hObject    handle to loadmaxprojections (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);hold off
cla reset;
handles.axes1.YTickLabel=[];
handles.axes1.XTickLabel=[];
box on
image.DirectoryProjections = { {'uigetdir(''.'')'} };
image.MicroscopePosition_Andor = 0;
image.channel = 0;
image.delta_t_in_minutes = 0;
%image.contrastLims=[0 0];
image = StructDlg(image);
handles.delta_t = image.delta_t_in_minutes;
handles.chan = image.channel;
handles.dir = image.DirectoryProjections;
handles.pos = image.MicroscopePosition_Andor;
%handles.contrastLims=image.contrastLims;% user will input initial contrast limits to show the image
handles.currT = 1;
%here load the correct channel max projections data
[multi_img,nt] = loadProjections(handles);
 handles.szt =nt;
 handles.allprojections =multi_img; 
% now display the images
imshow(multi_img(:,:,1),[0 3500],'Parent',handles.axes1);% todo: deal with contrast limits smarter
set(handles.slider1, 'Min', 1); 
set(handles.slider1, 'Max', nt);
set(handles.slider1, 'Value', 1); %init at timept 1
set(handles.slider1, 'SliderStep', [1 / (size(handles.allprojections,3)) ,...
    10 / (size(handles.allprojections,3)) ]);%
%here make the loadtracks button visible
handles.loadtrackingdata.Visible = 'On';
handles.counter= 0; % to increment the number of chosen cells to track
% if want to look at fluor data but no tracking data
handles.tp_notrackingdata = [];
handles.loadMatfile.Visible = 'On';
handles.qFluorescenceCheckBox.Visible = 'On';
handles.tracktoplot = [];% default track; gets reset onece the actual track is selected
handles.tracks = [];
handles.matfile = [];
handles.ClearButton.Visible = 'on';
handles.ResetImageSize.Visible = 'on';

guidata(hObject, handles);

% --- Executes on button press in loadtrackingdata.
function loadtrackingdata_Callback(hObject, eventdata, handles)
% hObject    handle to loadtrackingdata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.DisplayTotalTracksField,'String','N');
[tracks.Ilastikfile_TrackingData_h5, pathname] =uigetfile('*.h5','File Selector');%uigetfile('All Files (*.*)','File Selector');%{ {'uigetfile(''.'')'} }
handles.ShowTrackIDsState = 1; % initialize the state of the counter of how many times the button was pressed
handles.ShowFullTrack_State = 1;
tracks = StructDlg(tracks);
handles.ilastikfile = fullfile(pathname,tracks.Ilastikfile_TrackingData_h5);
[coordintime] = extracktIlastikTracks(handles);
handles.tracks = coordintime;
handles.colormap = prism; % default map
% these default values are only used if user clickes slider without
% selecting the trackid 
handles.tracktoplot = [];% default track; gets reset onece the actual track is selected
handles.matfile = [];% if no matfile,only want to see tracking data
handles.InputTrackID_toLocateOnImage.Visible = 'on';
handles.ResetImageSize.Visible = 'on';
set(handles.InputTrackID_toLocateOnImage,'String','LocateTrackID');
handles.total_tracks = size(coordintime,2);
handles.DisplayTotalTracksField.Visible = 'On';
handles.loadMatfile.Visible = 'On';
handles.qFluorescenceCheckBox.Visible = 'On';
handles.ShowTrackIDs.Visible = 'On';
set(handles.ShowTrackIDs,'String','ShowTrackIDs');
set(handles.DisplayTotalTracksField,'String',num2str(handles.total_tracks));
handles.TotalTracksStaticText.Visible = 'On';
handles.ClearButton.Visible = 'On';
% handles.counter= 0; % to increment the number of chosen cells to track
% show the first frame of the movie
UpdateImage(handles,1);hold on
disp('DONE');
guidata(hObject, handles);

function [multi_img,nt] = loadProjections(handles)
ff1 = readAndorDirectory(handles.dir); 
nucmoviefile = getAndorFileName(ff1,ff1.p(handles.pos),[],[],ff1.w(handles.chan));%getAndorFileName(files,pos,time,z,w)
nreader = bfGetReader(nucmoviefile);
nt = nreader.getSizeT;% get the number of frames
disp(['OPENED ' nucmoviefile ]);
proj = bfopen(nucmoviefile);
multi_img = zeros(size(proj{1}{1},1),size(proj{1}{1},2),nt);% initialize to image size and the second dim is the number of time points
for jj=1:nt
multi_img(:,:,jj) = (proj{1}{jj});
end

function UpdateImage(handles,time)
displayedImg=imshow(handles.allprojections(:,:,time),[500 3000],'Parent',handles.axes1);%hold on
set(displayedImg, 'ButtonDownFcn', @(src, evnt)IDclickedcell(src, evnt, handles));% @(src, evnt)IDclickedcell(src, evnt, handles)

function [coordintime] = extracktIlastikTracks(handles)
% here fully use Ilastik Automatic tracking 
clear coordintime
handles.pos = 1;%1
handles.chan = 1;% ff.w(chan): 0 - cfp cells; 1- nuc marker of other cell type
[coordintime] = process_ilastik_trackAN(handles.ilastikfile);

function DisplayTotalTracksField_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayTotalTracksField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DisplayTotalTracksField as text
%        str2double(get(hObject,'String')) returns contents of DisplayTotalTracksField as a double


% --- Executes during object creation, after setting all properties.
function DisplayTotalTracksField_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DisplayTotalTracksField (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plotpressed = 0;
handles.counter = 0;
handles.clickedX_notrackinfo = [];
handles.clickedY_notrackinfo = [];
handles.RandomcellID = [];
handles.quantifyselected = [];
handles.displayedFluorData = [];
handles.fulltrack = [];
handles.plotfluorinput = [];
handles.ShowTrackIDsState = 1;
handles.ShowFullTrack_State = 1;
handles.show_qFluorData = 0;
handles.validtrack = [];
handles.tracktoplot = [];
handles.trackcolor = [];
handles.clickedcellID=[];
handles.clickedcellX=[];
handles.clickedcellY=[];
handles.tp_notrackingdata = [];

% when there's no tracking data, the qFluor stays checked
if isempty(handles.tracks)
  set(handles.qFluorescenceCheckBox,'Value',1);
  handles.show_qFluorData = 1;
  handles.quantifyselected = 1;
else
set(handles.qFluorescenceCheckBox,'Value',0);
handles.show_qFluorData_State = 1;   
end

handles.SaveTrackIDs.Visible = 'Off';
handles.ResetImageSize.Visible = 'on';
handles.PlotFluorData_SelectedCells.Visible = 'Off';
handles.YlabelString.Visible = 'off';
set(handles.YlabelString,'String','PlotFluor Ylabel');
handles.InputTrackID_toLocateOnImage.Visible = 'on';
set(handles.InputTrackID_toLocateOnImage,'String','LocateTrackID');
handles.ShowTrackIDs.Visible = 'off';
set(handles.ShowTrackIDs,'String','ShowTrackIDs');
handles.ShowFullTrack.Visible = 'Off';
axes(handles.axes1);hold off
cla reset;
handles.axes1.YTickLabel=[];
handles.axes1.XTickLabel=[];
box on
delete(handles.axes1.Children)
set(handles.slider1, 'Value', 1); 
UpdateImage(handles,1);% show first frame when clear button is pressed
title(['Current time point  ' num2str(1*handles.delta_t/60)   'hrs since start; Frame(' num2str(1) ')']);hold on
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function TotalTracksStaticText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TotalTracksStaticText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in ShowFullTrack.
function ShowFullTrack_Callback(hObject, eventdata, handles)
% hObject    handle to ShowFullTrack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%handles.tracktoplot = str2double(get(handles.ReadTrackID, 'String'));
if isempty(handles.tracktoplot)
    disp('No tracks chosen to display. Click on a cell to choose')
    return
end
handles.ShowFullTrack_State = handles.ShowFullTrack_State+1;
xyt = handles.tracks;
handles.currT = int32(get(handles.slider1, 'Value'));
%disp(handles.ShowFullTrack_State);
if mod(handles.ShowFullTrack_State,2) >0 % if odd
set(handles.ShowFullTrack,'String','ShowFullTrack');
% keep only the cell label at current slider position
toview = size(handles.tracktoplot,2);
xyt_curr = struct;
for jj=1:toview
   xyt_curr(jj).fulltrack = xyt(handles.tracktoplot(jj)).dat;  
end
handles.fulltrack = xyt_curr; 
UpdateImage(handles,handles.currT);hold on
colormap = handles.colormap;
for k=1:toview
 firstT(k)= min(nonzeros(xyt_curr(k).fulltrack(:,3)));
 lastT(k) =  max(nonzeros(xyt_curr(k).fulltrack(:,3)));
if  (firstT(k)<=handles.currT) && (handles.currT<=lastT(k))% if this track has an XY at this time point
hplot2 = plot(xyt_curr(k).fulltrack(handles.currT,1),xyt_curr(k).fulltrack(handles.currT,2)...
    ,'p','MarkerFaceColor',colormap(handles.trackcolor(k),:),'MarkerEdgeColor',colormap(handles.trackcolor(k),:),'MarkerSize',5,'LineWidth',1);hold on
text(xyt_curr(k).fulltrack(handles.currT,1)+5,xyt_curr(k).fulltrack(handles.currT,2)+5,num2str(handles.tracktoplot(k)),'Color',colormap(handles.trackcolor(k),:),'FontSize',12);hold on
end
end
%title(['Frame ' num2str(handles.currT) ]);
end
if mod(handles.ShowFullTrack_State,2) ==0 % if even
set(handles.ShowFullTrack,'String','HideFullTrack');
toview = size(handles.tracktoplot,2);
xyt_curr = struct;
for jj=1:toview
   xyt_curr(jj).fulltrack = xyt(handles.tracktoplot(jj)).dat;  
end
handles.fulltrack = xyt_curr; 
handles.currT = int32(get(handles.slider1, 'Value'));
UpdateImage(handles,handles.currT);hold on
colormap = handles.colormap;
for k=1:toview
 firstT(k)= min(nonzeros(xyt_curr(k).fulltrack(:,3))); 
 lastT(k) =  max(nonzeros(xyt_curr(k).fulltrack(:,3)));
if  (firstT(k)<=handles.currT) && (handles.currT<=lastT(k)) % if this track has an XY at this time point
hplot2 = plot(xyt_curr(k).fulltrack(firstT(k)+1:end-1,1),xyt_curr(k).fulltrack(firstT(k)+1:end-1,2)...
    ,'-*','Color',colormap(handles.trackcolor(k),:),'MarkerSize',3,'LineWidth',1);hold on
plot(xyt_curr(k).fulltrack(firstT(k),1),xyt_curr(k).fulltrack(firstT(k),2)...
    ,'p','MarkerFaceColor','r','MarkerEdgeColor','w','MarkerSize',6,'LineWidth',1);hold on
text(xyt_curr(k).fulltrack(firstT(k),1)+8,xyt_curr(k).fulltrack(firstT(k),2)+5,'Start','Color','r','FontSize',8);hold on
plot(xyt_curr(k).fulltrack(end,1),xyt_curr(k).fulltrack(end,2)...
    ,'p','MarkerFaceColor','g','MarkerEdgeColor','w','MarkerSize',6,'LineWidth',1);hold on
text(xyt_curr(k).fulltrack(end,1)+8,xyt_curr(k).fulltrack(end,2)+5,'End','Color','g','FontSize',8);hold on
end
end
%tp = min(firstT);
%title('Showing Frame 1 ');
%set(handles.slider1, 'Value', 1); 
end
guidata(hObject, handles);


% --- Executes on button press in SaveTrackIDs.
function SaveTrackIDs_Callback(hObject, eventdata, handles)
% hObject    handle to SaveTrackIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
platf = ispc;% see if PC or MAC to use the '/' correctly for the save path
goodtracks = handles.validtrack;
if ~isempty(handles.tracks)
if platf == 1
    fname = [num2str(handles.dir) '\SelectedTrackIDsData_f' num2str(handles.pos-1) '_chan_w000' num2str(handles.chan-1) '.mat' ];
else
    fname = [num2str(handles.dir) '/SelectedTrackIDsData_f' num2str(handles.pos-1) '_chan_w000' num2str(handles.chan-1) '.mat' ];
end
% only the clicked cells are quntified and only
% they will be saved (no manual input of trackIDs to save)
if ~isempty(handles.matfile) && (handles.show_qFluorData==1)
    fluor_selectedtracks = struct;
    for jj=1:size(handles.quantifyselected,2)
    fluor_selectedtracks(jj).fluor = cat(1,handles.quantifyselected(jj).selected.dat.fluor); %handles.quantifyselected(jj).selected.dat
    fluor_selectedtracks(jj).xy = cat(1,handles.quantifyselected(jj).selected.dat.xy);
    fluor_selectedtracks(jj).trackID =handles.quantifyselected(jj).selected.dat(1).track;
    end
    save(fname,'goodtracks','fluor_selectedtracks');disp('Saved selected trackIDs and their fluorescence quantification to file');
else
    save(fname,'goodtracks');
    disp('Saved selected trackIDs to file');
end
else % if no tracking data was used    
    if platf == 1
    fname = [num2str(handles.dir) '\SelectedCellsData_randomTrackIDs_f' num2str(handles.pos-1) '_chan_w000' num2str(handles.chan-1) '.mat' ];
else
    fname = [num2str(handles.dir) '/SelectedTrackIDsData_randomTrackIDs_f' num2str(handles.pos-1) '_chan_w000' num2str(handles.chan-1) '.mat' ];
end
% only the clicked cells are quntified and only
% they will be saved (no manual input of trackIDs to save)
if ~isempty(handles.matfile) && (handles.show_qFluorData==1)
    fluor_selectedtracks = struct;
    for jj=1:size(handles.quantifyselected,2)
    fluor_selectedtracks(jj).fluor = cat(1,handles.quantifyselected(jj).selected.dat.fluor); %handles.quantifyselected(jj).selected.dat
    fluor_selectedtracks(jj).xy = cat(1,handles.quantifyselected(jj).selected.dat.xy);
    fluor_selectedtracks(jj).trackID =cat(1,handles.quantifyselected(jj).selected.dat.track);
    fluor_selectedtracks(jj).tp =cat(1,handles.quantifyselected(jj).selected.dat.tp);
    end
    save(fname,'fluor_selectedtracks');disp('Saved selected trackIDs and their fluorescence quantification to file');
end
end

guidata(hObject, handles);

% --- Executes on button press in ShowTrackIDs.
function ShowTrackIDs_Callback(hObject, eventdata, handles)
% hObject    handle to ShowTrackIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ShowTrackIDsState = handles.ShowTrackIDsState+1; % anytime the button is clicked, this variable is incremented
%get the slider position
handles.currT = int32(get(handles.slider1, 'Value'));
%UpdateImage(handles,1);hold on
%disp(handles.ShowTrackIDsState);
if mod(handles.ShowTrackIDsState,2) >0 % if odd
set(handles.ShowTrackIDs,'String','ShowTrackIDs');
UpdateImage(handles,handles.currT);hold on
end
if mod(handles.ShowTrackIDsState,2) == 0 % if even
set(handles.ShowTrackIDs,'String','HideTrackIDs');
% code to plot all the track IDs on the image
xyt = handles.tracks;
allIDdata = [];
for jj=1:size(xyt,2)
    firstT = min(nonzeros(xyt(jj).dat(:,3)));%does track exist at the currentT?
    lastT = max(nonzeros(xyt(jj).dat(:,3))); %does track exist at the currentT
    if (firstT<= handles.currT) &&  (handles.currT<=lastT)
        allIDdata(jj,1:2)  = xyt(jj).dat(handles.currT,1:2);
        allIDdata(jj,3)=jj;    % trackID label
    end
 end
hdata  = plot(allIDdata(:,1),allIDdata(:,2)...
    ,'*m','MarkerSize',3,'LineWidth',1);hold on
text(allIDdata(:,1)+3*ones(size(allIDdata(:,1))),allIDdata(:,2)+5*ones(size(allIDdata(:,1))),num2str(allIDdata(:,3)),'FontSize',7,'Color','g');
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ShowTrackIDs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ShowTrackIDs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.ShowFullTrack.Visible = 'On';


function IDclickedcell(src,evnt,handles)
handles = guidata(src);
handles.colormap = prism;     
clickedXY = get(handles.axes1, 'CurrentPoint');
xy = [clickedXY(1,1) clickedXY(1,2)];
handles.currT = int32(get(handles.slider1, 'Value'));

% store the xy for when no tracked data is used, only the fluor
% find the trackID of the clicked cell 

handles.counter = handles.counter+1;% any tme a new track is input, the counter is incremented    
if ~isempty(handles.tracks)
handles.InputTrackID_toLocateOnImage.Visible = 'on';
handles.ShowFullTrack.Visible = 'On';
handles.ShowTrackIDs.Visible = 'On';
set(handles.ShowFullTrack,'String','ShowFullTrack');
set(handles.ShowTrackIDs,'String','ShowTrackIDs');    
xyt = handles.tracks;
allIDdata = [];
handles.currT = int32(get(handles.slider1, 'Value'));
for jj=1:size(xyt,2)
    firstT = min(nonzeros(xyt(jj).dat(:,3)));%does track exist at the currentT?
    lastT = max(nonzeros(xyt(jj).dat(:,3))); %does track exist at the currentT
    if (firstT<= handles.currT) &&  (handles.currT<=lastT)
        allIDdata(jj,1:2)  = xyt(jj).dat(handles.currT,1:2);
        allIDdata(jj,3)=jj;    % trackID label
    end
end
closestcell= ipdm(xy,allIDdata(:,1:2), 'Subset', 'NearestNeighbor', 'Result', 'Structure');
% here put the readtrackID code and get rid of the button, such that
% theclick populates the totrack 
readinput= allIDdata(closestcell.columnindex,3);% trackID as read from the click
handles.tracktoplot(handles.counter) = readinput;
handles.colormap = prism;      
handles.PlotTrackID.Visible = 'On';
handles.ShowFullTrack.Visible = 'On';
colormap = handles.colormap;
randcolor = randi(size(colormap,1));%randi(size(colormap,1))     select the track label color once
handles.trackcolor(handles.counter) = randcolor;
handles.clickedcellID(handles.counter) = allIDdata(closestcell.columnindex,3);
handles.clickedcellX(handles.counter)=allIDdata(closestcell.columnindex,1);
handles.clickedcellY(handles.counter)=allIDdata(closestcell.columnindex,2);
% if the fist cell is uncklicked
if handles.counter ==2 && (handles.tracktoplot(handles.counter) == handles.tracktoplot(handles.counter-1))
  handles.counter = 0;
  handles.tracktoplot = [];
  handles.trackcolor = [];  
  handles.clickedcellID = [];
  handles.clickedcellY = [];
  handles.clickedcellX = [];
  handles.currT = int32(get(handles.slider1, 'Value'));
  UpdateImage(handles,handles.currT); 
 
end
% here if the same cell was clicekd again (consecutively), make the cell label dissapear; 
if (handles.counter > 2) && (handles.tracktoplot(handles.counter) == handles.tracktoplot(handles.counter-1))
  handles.counter = handles.counter-2;
  handles.tracktoplot = (handles.tracktoplot(1:end-2));  
  handles.trackcolor=handles.trackcolor(1:handles.counter);  
  handles.clickedcellID = handles.clickedcellID(1:(handles.counter));
  handles.clickedcellX = handles.clickedcellX(1:(handles.counter));
  handles.clickedcellY=  handles.clickedcellY(1:(handles.counter)); 
  handles.currT = int32(get(handles.slider1, 'Value'));
  UpdateImage(handles,handles.currT); 
  end
% if the cell was unclicked later on,after selecting some more cells
if (handles.counter > 2) && any(handles.tracktoplot(handles.counter) == (handles.tracktoplot(1:handles.counter-1)))
  [~,c]=find(handles.tracktoplot(handles.counter) ~= handles.tracktoplot);% any trackIDs that don't match the clicked one
 % if there's no tracking data
  handles.counter = handles.counter-2;
  handles.tracktoplot = handles.tracktoplot(c);  
  handles.trackcolor= handles.trackcolor(c);  
  handles.clickedcellID = handles.clickedcellID(c);
  handles.clickedcellX = handles.clickedcellX(c);
  handles.clickedcellY=  handles.clickedcellY(c);
  handles.currT = int32(get(handles.slider1, 'Value'));
  UpdateImage(handles,handles.currT);
  
end
else
handles.currT = int32(get(handles.slider1, 'Value'));
colormap = handles.colormap;
randcolor = randi(size(colormap,1));%randi(size(colormap,1))     select the track label color once
handles.trackcolor(handles.counter) = randcolor;
handles.tracktoplot(handles.counter) = randi(1000);%generate a random number to assign to the untrackd cell
handles.RandomcellID(handles.counter) = handles.tracktoplot(handles.counter);
handles.clickedX_notrackinfo(handles.counter) = round(xy(1));
handles.clickedY_notrackinfo(handles.counter) = round(xy(2));
disp(num2str(handles.tracktoplot(handles.counter)));
handles.tp_notrackingdata(handles.counter) = double(handles.currT);

end
%----------------------
% tdata.String = [];
handles.currT = int32(get(handles.slider1, 'Value'));
UpdateImage(handles,handles.currT);
if ~isempty(handles.tracks)
    fluor_selectedtracks = struct;
    handles.quantifyselected = [];
    for jj=1:handles.counter%size(handles.tracktoplot,2)        
        hdata  = plot(handles.clickedcellX(jj),handles.clickedcellY(jj)...
            ,'p','MarkerFaceColor',colormap(handles.trackcolor(jj),:),'MarkerEdgeColor',colormap(handles.trackcolor(jj),:),'MarkerSize',5,'LineWidth',1);hold on
        tdata = text(handles.clickedcellX(jj)+5,handles.clickedcellY(jj)+5,...
            num2str(handles.clickedcellID(jj)),'FontSize',12,'Color',colormap(handles.trackcolor(jj),:));

        if ~isempty(handles.matfile) && (handles.show_qFluorData == 1)
            % the function below will get the fluoresence quantification data for each
            % clicked cell for all time points availabe in the .mat file
            [fluor_dat]=getfluordata(handles,handles.tracktoplot(jj));%
            % need then to store it in handles and use it at slider motion if the .mat
            % file is loaded
            fluor_selectedtracks(jj).dat=fluor_dat;
            handles.quantifyselected(jj).selected= fluor_selectedtracks(jj);
            
        end
    end
     disp(['Cell tracks to watch simultaneously: ' num2str(handles.counter) ]);

else
     % only plot markers for the slider time point
        [~,c4]=find(handles.tp_notrackingdata == double(handles.currT));        
    for jj=1:size(c4,2)        
            hdata  = plot(handles.clickedX_notrackinfo(c4(jj)),handles.clickedY_notrackinfo(c4(jj))...
            ,'p','MarkerFaceColor',colormap(handles.trackcolor(c4(jj)),:),'MarkerEdgeColor',colormap(handles.trackcolor(c4(jj)),:),'MarkerSize',5,'LineWidth',1);hold on
                tdata = text(handles.clickedX_notrackinfo(c4(jj))+5,handles.clickedY_notrackinfo(c4(jj))+5,...
            num2str(handles.RandomcellID(c4(jj))),'FontSize',12,'Color',colormap(handles.trackcolor(c4(jj)),:));hold on
       % this will find the best match to the clicked centroid in the mat
       % file
      
        
    end
    [xy_centroid]=find_cellinmatfile(handles.clickedX_notrackinfo(end),handles.clickedY_notrackinfo(end),handles);
%             % plot the coordinates of the actual cell that was matched to
%             % clicked cell in the fluor matfile
            hdata  = plot(round(xy_centroid(1)),round(xy_centroid(2))...
            ,'p','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',5,'LineWidth',1);hold on

    
 disp(['Total cells selected: ' num2str(handles.counter) ]);
  
end
if ~isempty(handles.matfile) && (handles.show_qFluorData == 1)
    handles.PlotFluorData_SelectedCells.Visible = 'On';
     if isempty(handles.displayedFluorData)
    handles.YlabelString.Visible = 'on';
    set(handles.YlabelString,'String','PlotFluor Y label');
    end
end 

%handles.validtrack = handles.clickedcellID;
handles.validtrack = handles.tracktoplot;
handles.SaveTrackIDs.Visible = 'On';
guidata(src, handles);

% --- Executes on button press in loadMatfile.
function loadMatfile_Callback(hObject, eventdata, handles)
% hObject    handle to loadMatfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% here will need to have some checks on whether the file has the info on
% each time point quntified or only the last time point
handles.show_qFluorData = 0;
handles.show_qFluorData_State =1;
[tracks.fluorquantified, pathname_mat] =uigetfile('*.mat','File Selector');% 
tracks = StructDlg(tracks);
handles.matfile = fullfile(pathname_mat,tracks.fluorquantified);% the output .mat file of the Idse's analysis code, with the 'positions' containing all the data
handles.quantifyselected = [];
fluor_selectedtracks = struct;
if isempty(handles.tracks) 
    handles.plotpressed = 0;
set(handles.qFluorescenceCheckBox,'Value',1);
handles.show_qFluorData = 1;
handles.displayedFluorData = [];
handles.quantifyselected = 1;
end
disp('Done');
guidata(hObject, handles);

function [closestcell]=getfluordata(handles,readinput)
% helper function, called within slider callback
load(handles.matfile);

if ~isempty(handles.tracks)
xyt = handles.tracks; 
xyt_curr = struct;
if ~isempty(handles.tracktoplot)
   xyt_curr.fulltrack = xyt(readinput).dat;   
end
jj = handles.pos;
closestcell_tmp=[];
closestcell = struct;
tpt = size(xyt_curr.fulltrack,1);% how mny time points were tracked
bg = positions(handles.pos).timeTraces.background;%
 tracktoquantify = readinput;
 matched = struct;
 tomatch = struct;
     for ii=1:tpt% time points; TODO make sure the ii and the third column of the tracked data are the same time point
        if ~isempty(positions(jj).cellData)
            tomatch(ii).xy = xyt_curr.fulltrack(ii,1:2);% TODO: make sure this time point data exists
            % find the tracked cell among the fluorQ data (it has no
            % trackID, so need to find the closest cell to the tracked cell at each time point)
            allcells_givenT = positions(jj).cellData(ii).XY;
            closestcell_tmp= ipdm(tomatch(ii).xy,allcells_givenT, 'Subset', 'NearestNeighbor', 'Result', 'Structure');
            % get fluor data for all cells at that time point
            rmbg = ones(size(positions(jj).cellData(ii).nucLevel))*bg(ii);% get bg vector; it is different for each time point
            tmp = (positions(jj).cellData(ii).nucLevel(:,1)-rmbg)./(positions(jj).cellData(ii).cytLevel(:,1)-rmbg); % subtrct bg from nuc and cyto green levels
            closestcell(ii).fluor = tmp(closestcell_tmp.columnindex,:);
            closestcell(ii).xy = allcells_givenT(closestcell_tmp.columnindex,:);% these are the coordinates from the .mat file, not tracking
            closestcell(ii).distance = closestcell_tmp.distance;
            closestcell(ii).track = readinput;
        end
     end
else
     % if want to look at q fluor in untracked data
jj = handles.pos;
closestcell_tmp=[];
closestcell = struct;
tpts = handles.tp_notrackingdata;% 
bg = positions(handles.pos).timeTraces.background;%
 tracktoquantify = readinput;
 matched = struct;
 tomatch = struct;
     for ii=1:size(tpts,2)
         % current time point;
        if ~isempty(positions(jj).cellData) 
            %[~,c3]=find(handles.RandomcellID == readinput);
            tomatch(ii).xy = cat(2,handles.clickedX_notrackinfo(ii),handles.clickedY_notrackinfo(ii));
            % trackID, so need to find the closest cell to the clicked cell at given time point)
            allcells_givenT = positions(jj).cellData(tpts(ii)).XY;
            closestcell_tmp= ipdm(tomatch(ii).xy,allcells_givenT, 'Subset', 'NearestNeighbor', 'Result', 'Structure');
            % get fluor data for all cells at that time point
            rmbg = ones(size(positions(jj).cellData(tpts(ii)).nucLevel))*bg(tpts(ii));% get bg vector; it is different for each time point
            tmp = (positions(jj).cellData(tpts(ii)).nucLevel(:,1)-rmbg)./(positions(jj).cellData(tpts(ii)).cytLevel(:,1)-rmbg); % subtrct bg from nuc and cyto green levels
            closestcell(ii).fluor = tmp(closestcell_tmp.columnindex,:);
            closestcell(ii).xy = allcells_givenT(closestcell_tmp.columnindex,:);% these are the coordinates from the .mat file, not tracking
            closestcell(ii).distance = closestcell_tmp.distance;
            closestcell(ii).track = handles.RandomcellID(ii);
            closestcell(ii).tp = tpts(ii);
           % disp(['Centroid found in qFluor data at d =  :' num2str(closestcell_tmp.distance) 'pxls from clicked spot' ])% find the clicked cell among the fluorQ data (it has no
            
        end
     end 
    end

% --- Executes on button press in qFluorescenceCheckBox.
function qFluorescenceCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to qFluorescenceCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.matfile)
    disp('Load Mat File first');
    set(handles.qFluorescenceCheckBox,'Value',0);
    return
end
 handles.show_qFluorData_State= handles.show_qFluorData_State+1;
if mod(handles.show_qFluorData_State,2)==0
    handles.show_qFluorData = 1;
end
if mod(handles.show_qFluorData_State,2)>0
    handles.show_qFluorData  = 0;
end
%disp(handles.show_qFluorData);  
    handles.displayedFluorData = [];
guidata(hObject, handles);

% --- Executes on button press in PlotFluorData_SelectedCells.
function PlotFluorData_SelectedCells_Callback(hObject, eventdata, handles)
% hObject    handle to PlotFluorData_SelectedCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fluorallt = struct;
if isempty(handles.displayedFluorData)
    disp('Input the ylabel for the fluorescence data plot')
    return
else
colormap = handles.colormap;
if isempty(handles.quantifyselected)
    disp('Select cells after clicking the checkbox')
    return
end
%TODO: add the different symbols, since sometimes the tracks have the same
%color
legendstr = cell(size(handles.quantifyselected,2),1);
symbolstring = {'x','p','d','o','s','*','x','p','d','o','s','*'};%
if ~isempty(handles.tracks)
if size(symbolstring,2)<size(handles.quantifyselected,2)
    disp('Plotting all data using the same marker')
    for x=1:size(handles.quantifyselected,2)        
        fluorallt(x).intensity= cat(1,handles.quantifyselected(x).selected.dat.fluor);
        fluorallt(x).coord = round(cat(1,handles.quantifyselected(x).selected.dat.xy));
        figure(1),plot(fluorallt(x).intensity,'Color',...
            colormap(handles.trackcolor(x),:),'LineWidth',1);hold on; box on
        limitsX(x) = size(fluorallt(x).intensity,1);
        limitsY(x) = max(fluorallt(x).intensity);
        legendstr{x} = num2str(handles.tracktoplot(x));
    end
else
for x=1:size(handles.quantifyselected,2)    
 fluorallt(x).intensity= cat(1,handles.quantifyselected(x).selected.dat.fluor);
 fluorallt(x).coord = round(cat(1,handles.quantifyselected(x).selected.dat.xy));
 figure(1),plot(fluorallt(x).intensity,'Color',...
     colormap(handles.trackcolor(x),:),'Marker',symbolstring{x},'LineWidth',2);hold on; box on 
 limitsX(x) = size(fluorallt(x).intensity,1);
 limitsY(x) = max(fluorallt(x).intensity);
 legendstr{x} = num2str(handles.tracktoplot(x));
end
end
else
handles.plotpressed=1;
fluor_selectedtracks = struct;
handles.quantifyselected = [];  
if ~isempty(handles.matfile) && (handles.show_qFluorData == 1)
    % the function below will get the fluoresence quantification data for each
    % clicked cell for all time points availabe in the .mat file
    [fluor_dat]=getfluordata(handles,[]);  %jj
    % need then to store it in handles and use it at slider motion if the .mat
    % file is loaded
    for jj=1
    fluor_selectedtracks(jj).dat=fluor_dat;%jj
    handles.quantifyselected(jj).selected= fluor_selectedtracks(jj);%jj
    end
end

    handles.currT = int32(get(handles.slider1, 'Value'));
    if size(symbolstring,2)<size(handles.quantifyselected(1).selected.dat,2)
    disp('Plotting all data using the same marker')
 fluorallt.intensity= cat(1,handles.quantifyselected(1).selected.dat.fluor);
 fluorallt.coord = round(cat(1,handles.quantifyselected(1).selected.dat.xy));
 fluorallt.selectedtimes = round(cat(1,handles.quantifyselected(1).selected.dat.tp));
 newsz=size(nonzeros(isfinite(fluorallt.intensity)),1);
 new_t = fluorallt.selectedtimes(isfinite(fluorallt.intensity));
 new_I =fluorallt.intensity(isfinite(fluorallt.intensity));
 new_c =colormap(handles.trackcolor(isfinite(fluorallt.intensity)),:);
 legendstr = cell(1,newsz)  ;
  for w=1:newsz
 figure(1),plot(new_t(w),new_I(w),'Color',...
     new_c(w,:),'Marker','p','LineWidth',2);hold on; box on 
  end
 limitsX = size(new_I,1);
 limitsY = max(new_I); 
else

 fluorallt.intensity= cat(1,handles.quantifyselected(1).selected.dat.fluor);
 fluorallt.coord = round(cat(1,handles.quantifyselected(1).selected.dat.xy));
 fluorallt.selectedtimes = round(cat(1,handles.quantifyselected(1).selected.dat.tp));
 newsz=size(nonzeros(isfinite(fluorallt.intensity)),1);
 new_t = fluorallt.selectedtimes(isfinite(fluorallt.intensity));
 new_I =fluorallt.intensity(isfinite(fluorallt.intensity));
 new_c =colormap(handles.trackcolor(isfinite(fluorallt.intensity)),:);
 legendstr = cell(1,newsz)  ;
  for w=1:newsz
 figure(1),plot(new_t(w),new_I(w),'Color',...
     new_c(w,:),'Marker',symbolstring{w},'LineWidth',2);hold on; box on 
  end
 % figure(1), errorbar()
 limitsX = size(new_I,1);
 limitsY = max(new_I); 
 
end
end
f = figure(1);
if ~isempty(handles.tracks)
title('Selected cells only, Legend: trackIDs from Ilastik');
legend(legendstr,'Location','NorthWest');
if isfinite(max(limitsY)) && isfinite(max(limitsX)) 
ylim([0 (max(limitsY)+0.5)]);
xlim([0 (max(limitsX)+1)]);
f.CurrentAxes.XTick = 1:9:max(limitsX);
f.CurrentAxes.XTickLabel = (1:9:max(limitsX))*handles.delta_t/60;
else
    disp('One of the data points is NaN')    
end
else
  ylim([0 (max(limitsY)+0.5)]);
  xlim([0 handles.szt+1]);
  f.CurrentAxes.XTick = 1:17:handles.szt;
  f.CurrentAxes.XTickLabel = (1:17:handles.szt)*handles.delta_t/60;  
  title('Selected cells only, no tracking data ');
end
f.CurrentAxes.FontSize = 12;
f.CurrentAxes.LineWidth = 2;
xlabel('Time, hrs');    
ylabel(handles.displayedFluorData);
end
guidata(hObject, handles);


function InputTrackID_toLocateOnImage_Callback(hObject, eventdata, handles)
% hObject    handle to InputTrackID_toLocateOnImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputTrackID_toLocateOnImage as text
%        str2double(get(hObject,'String')) returns contents of InputTrackID_toLocateOnImage as a double
handles.ShowFullTrack.Visible = 'On';
set(handles.ShowFullTrack,'String','ShowFullTrack');
set(handles.ShowTrackIDs,'String','ShowTrackIDs');
handles.SaveTrackIDs.Visible = 'On';
handles.counter = handles.counter+1;% any tme a new track is input, the counter is incremented    
inputID = str2double(get(handles.InputTrackID_toLocateOnImage, 'String'));
handles.tracktoplot(handles.counter) = inputID;
if (inputID > handles.total_tracks)
    disp('This number exceeds total number of tracks');
    set(handles.InputTrackID_toLocateOnImage,'String','');
    return
end
% here need to put some way of unclicking the manually input cell (not really useful, but may happen)

xyt = handles.tracks;    
handles.currT = int32(get(handles.slider1, 'Value'));
%set(handles.ShowTrackIDs,'String','ShowTrackIDs');
toview = size(handles.tracktoplot,2);
hold on;UpdateImage(handles,handles.currT);
if ~isempty(handles.tracktoplot)
xyt_curr = struct;sz_tmp = [];
for jj=1:toview
   xyt_curr(jj).fulltrack = xyt(handles.tracktoplot(jj)).dat;
   sz_tmp(jj) = size(xyt(handles.tracktoplot(jj)).dat,1);   
end
handles.fulltrack = xyt_curr;
colormap = handles.colormap;
randcolor =randi(size(colormap,1)) ;%randi(size(colormap,1))  select the track label color once here
handles.trackcolor(handles.counter) = randcolor;
for k =1:toview
 firstT(k)= min(nonzeros(xyt_curr(k).fulltrack(:,3)));
    if  handles.currT<=size(xyt_curr(k).fulltrack,1)% not all cells are tracked all the way
        if xyt_curr(k).fulltrack(handles.currT,1)>0 % tracks that were assigned after first time point, have zeros in the beginning so don't plot them
            plot(xyt_curr(k).fulltrack(handles.currT,1),xyt_curr(k).fulltrack(handles.currT,2),'p',...
                'MarkerFaceColor',colormap(handles.trackcolor(k),:),'MarkerEdgeColor',colormap(handles.trackcolor(k),:),'Markersize',2);hold on  
            text(handles.axes1,xyt_curr(k).fulltrack(handles.currT,1)+5,xyt_curr(k).fulltrack(handles.currT,2)+5,...
                num2str(handles.tracktoplot(k)),'Color',colormap(handles.trackcolor(k),:),'FontSize',12);hold on
        else
            disp(['Track ID ' num2str(handles.tracktoplot(k)) 'was picked up in frame' num2str(firstT(k))]);
        end

end
end
title(['Current time point  ' num2str((double(handles.currT)*handles.delta_t)/60)  ...
    'hrs since start; Frame(' num2str(handles.currT) ')']);hold on
end

if ~isempty(handles.matfile) && (handles.show_qFluorData == 1)
 handles.PlotFluorData_SelectedCells.Visible = 'On';
 if isempty(handles.displayedFluorData)
    handles.YlabelString.Visible = 'on';
    set(handles.YlabelString,'String','PlotFluor Y label');
    end

[fluor_dat]=getfluordata(handles,inputID);
fluor_selectedtracks(handles.counter).dat=fluor_dat;
handles.quantifyselected(handles.counter).selected= fluor_selectedtracks(handles.counter);
end
disp(['Cell tracks to watch simultaneously: ' num2str(handles.counter) ]);
set(handles.InputTrackID_toLocateOnImage,'String','');
handles.validtrack = handles.tracktoplot;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function InputTrackID_toLocateOnImage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputTrackID_toLocateOnImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function YlabelString_Callback(hObject, eventdata, handles)
% hObject    handle to YlabelString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of YlabelString as text
%        str2double(get(hObject,'String')) returns contents of YlabelString as a double
displayedData = (get(handles.YlabelString, 'String'));

handles.displayedFluorData = displayedData;
handles.YlabelString.Visible = 'off';
%set(handles.YlabelString,'String','');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function YlabelString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YlabelString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ResetImageSize.
function ResetImageSize_Callback(hObject, eventdata, handles)
% hObject    handle to ResetImageSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.axes1.XLim = [0.5000 1.0245e+03];
handles.axes1.YLim =  [0.5000 1.0245e+03];
guidata(hObject, handles);

function [xy_centroid]=find_cellinmatfile(x,y,handles)
load(handles.matfile);
jj = handles.pos;
closestcell_tmp=[];
closestcell = struct;
tpts = double(handles.currT);% 
tomatch = struct;
     for ii=tpts% current time point;         
        if ~isempty(positions(jj).cellData) 
            tomatch.xy = cat(2,x,y);
            allcells_givenT = positions(jj).cellData(ii).XY;
            closestcell_tmp= ipdm(tomatch.xy,allcells_givenT, 'Subset', 'NearestNeighbor', 'Result', 'Structure');
            closestcell.xy = allcells_givenT(closestcell_tmp.columnindex,:);% these are the coordinates from the .mat file, not tracking
            closestcell.distance = closestcell_tmp.distance;
          %  disp(['Centroid found in qFluor data at d =  :' num2str(closestcell_tmp.distance) 'pxls from clicked spot' ])% find the clicked cell among the fluorQ data (it has no
     end
     end 
     xy_centroid = closestcell.xy;

