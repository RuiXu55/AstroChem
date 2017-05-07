function varargout = GUI_sample(varargin)
% GUI_SAMPLE MATLAB code for GUI_sample.fig
%      GUI_SAMPLE, by itself, creates a new GUI_SAMPLE or raises the existing
%      singleton*.
%
%      H = GUI_SAMPLE returns the handle to a new GUI_SAMPLE or the handle to
%      the existing singleton*.
%
%      GUI_SAMPLE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SAMPLE.M with the given input arguments.
%
%      GUI_SAMPLE('Property','Value',...) creates a new GUI_SAMPLE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_sample_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_sample_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_sample

% Last Modified by GUIDE v2.5 12-Nov-2014 10:14:22

% Written by Rui Xu. Oct. 2014

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_sample_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_sample_OutputFcn, ...
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




% --- Executes just before GUI_sample is made visible.
function GUI_sample_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_sample (see VARARGIN)


%[sigma, species_name] = readfile(500);  
%%%%%%%%%%%%%%%%Read file%%%%%%%%%%%%%%%%%%

% Choose default command line output for GUI_sample
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_sample_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit_species_Callback(hObject, eventdata, handles)
% hObject    handle to edit_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_species as text
%        str2double(get(hObject,'String')) returns contents of edit_species as a double
if get(handles.add_complex,'value')==1
    complex_label =1;
else
    complex_label =0;
end

str = get(hObject,'string');
str_num1 = find_name(handles.species_name1,str);
str_num2 = find_name(handles.species_name2,str);
if str_num1 == -1 || str_num2 ==-1
    s = sprintf('No such kind of species in current chemistry network!');
    h = msgbox(s,'Warning Message','warn');
    return;
else
    axes(handles.axes1)
    semilogy(handles.z,handles.sigma1(:,str_num1)./handles.abundance_H1(:),'r-','LineWidth',2.5);
    temp_legend(1) = handles.species_name1(str_num1);
    hold on;
    if complex_label ==1
        semilogy(handles.z,handles.sigma4(:,str_num2)./handles.abundance_H4(:),'r--','LineWidth',2.5);
        temp_legend(2) = strcat(handles.species_name2(str_num2),'-cplx');
    end
    legend(temp_legend,'Location','Southeast');
    xlabel('Scale Height');
    ylabel('Relative Abundance /[H]');
    set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
    hold off;
    
    axes(handles.axes2)
    semilogy(handles.z,handles.sigma2(:,str_num1)./handles.abundance_H2(:),'r-','LineWidth',2.5);
    hold on;
    if complex_label ==1
        semilogy(handles.z,handles.sigma5(:,str_num2)./handles.abundance_H5(:),'r--','LineWidth',2.5);
    end    
    legend(temp_legend,'Location','Southeast');
    xlabel('Scale Height');
    ylabel('Relative Abundance /[H]');
    set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
    hold off;
    
    axes(handles.axes3)
    semilogy(handles.z,handles.sigma3(:,str_num1)./handles.abundance_H3(:),'r-','LineWidth',2.5);
    hold on;
    if complex_label ==1
        semilogy(handles.z,handles.sigma6(:,str_num2)./handles.abundance_H6(:),'r--','LineWidth',2.5);
    end
     legend(temp_legend,'Location','Southeast');
    xlabel('Scale Height');
    ylabel('Relative Abundance /[H]');
    set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
    hold off;
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_species_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Clear.
function Clear_Callback(hObject, eventdata, handles)
% hObject    handle to Clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
legend('off');
cla(handles.axes1);
axes(handles.axes2);
legend('off');
cla(handles.axes2);
axes(handles.axes3);
legend('off');
cla(handles.axes3);

% --- Executes on button press in Quit.
function Quit_Callback(hObject, eventdata, handles)
% hObject    handle to Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%close(gco);
delete(handles.figure1);

% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
legend('off');
cla(handles.axes1);
axes(handles.axes2);
legend('off');
cla(handles.axes2);
axes(handles.axes3);
legend('off');
cla(handles.axes3);

set(handles.edit_species,'string',' ');
set(handles.vertical_range,'string',' ');
set(handles.NetLabel,'string',' ');
set(handles.Scale_height,'string',' ');
set(handles.Vertical_num,'string',' ');
set(handles.H2O,'value',0);
set(handles.CO2,'value',0);
set(handles.CO,'value',0);
set(handles.OH,'value',0);
set(handles.H3O,'value',0);
set(handles.add_complex,'value',0);



% --- Executes on button press in Compare.
function [legendObject] = Compare_Callback(hObject, eventdata, handles)
% hObject    handle to Compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%num =handles.species_num;
if get(handles.add_complex,'value')==1
    complex_label =1;
else
    complex_label =0;
end

prompt = {'Name of the species'};
dlg_title = 'Input';
num_lines =1;
%def = {'20','hsv'};   %default value for the input
answer = inputdlg(prompt,dlg_title,num_lines);
% name is the compared species name
name = regexp(answer{1},'\s+','split');
name_label1 = find_multi_name(handles.species_name1, name);
name_label2 = find_multi_name(handles.species_name2, name);
num =length(name);
for i=1:num
    if name_label1(i) == -1 || name_label2(i) ==-1
        s = sprintf('No such kind of species in current chemistry network!');
        h = msgbox(s,'Warning Message','warn');
        return;
    end
end

axes(handles.axes1)
% set(0,'DefaultAxesColorOrder',[1 0 0;0 1 0; 0 0 1],...
%     'DefaultAxesLineStyleOrder','-|-.|--|:');
k =1;
map_inter =int16(fix(256/num));
for i=1:num
    semilogy(handles.z,handles.sigma1(:,name_label1(i))./handles.abundance_H1(:), ...
        'color',handles.mymap(map_inter*i,:),'LineStyle','-','LineWidth',2.5);
    temp_legend(k) = handles.species_name1(name_label1(i));
    k =k +1;
    hold all;
    if complex_label ==1
        semilogy(handles.z,handles.sigma4(:,name_label2(i))./handles.abundance_H4(:), ...
            'color',handles.mymap(map_inter*i,:),'LineStyle','--','LineWidth',2.5);
        temp_legend(k) = strcat(handles.species_name2(name_label2(i)),'-cplx');
        k = k+1;
    end
end
legendObject = legend(temp_legend,'Location','Southeast');
xlabel('Scale Height');
ylabel('Relative Abundance /[H]');
set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
% set(0,'DefaultAxesLineStyleOrder','remove');
% set(0,'DefaultAxesColorOrder','remove');
hold off;

axes(handles.axes2)

for i=1:num
    semilogy(handles.z,handles.sigma2(:,name_label1(i))./handles.abundance_H2(:), ...
        'color',handles.mymap(map_inter*i,:),'LineStyle','-','LineWidth',2.5);
    hold all;
    if complex_label ==1
        semilogy(handles.z,handles.sigma5(:,name_label2(i))./handles.abundance_H5(:), ...
            'color',handles.mymap(map_inter*i,:),'LineStyle','--','LineWidth',2.5);
    end
end
legend(temp_legend,'Location','Southeast');
xlabel('Scale Height');
set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
hold off;

axes(handles.axes3)

for i=1:num
    semilogy(handles.z,handles.sigma3(:,name_label1(i))./handles.abundance_H3(:), ...
        'color',handles.mymap(map_inter*i,:),'LineStyle','-','LineWidth',2.5);
    hold all;
    if complex_label ==1
        semilogy(handles.z,handles.sigma6(:,name_label2(i))./handles.abundance_H6(:), ...
        'color',handles.mymap(map_inter*i,:),'LineStyle','--','LineWidth',2.5);
    end
end
legend(temp_legend,'Location','Southeast');
xlabel('Scale Height');
set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
hold off;

set(legendObject,'Interpreter','none');
handles.legendObject = legendObject;
guidata(hObject,handles);

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over Compare.
function Compare_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes during object creation, after setting all properties.
function number_species_CreateFcn(hObject, eventdata, handles)
% hObject    handle to number_species (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in H2O.
function H2O_Callback(hObject, eventdata, handles)
% hObject    handle to H2O (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of H2O


% --- Executes on button press in CO2.
function CO2_Callback(hObject, eventdata, handles)
% hObject    handle to CO2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CO2


% --- Executes on button press in OH.
function OH_Callback(hObject, eventdata, handles)
% hObject    handle to OH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of OH


% --- Executes on button press in H3O.
function H3O_Callback(hObject, eventdata, handles)
% hObject    handle to H3O (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of H3O



% --- Executes on button press in plot.
function [legendObject] = plot_Callback(hObject, eventdata, handles)
% hObject    handle to plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%num indicate whether a certain species is selected
if get(handles.H2O,'value')==1
    num(1) = find_name(handles.species_name1,'H2O');
    num(6) = find_name(handles.species_name2,'H2O');
else
    num(1) = -0.5;         %num =-0.5 shows that do not plot this species         
end
if get(handles.CO2,'value')==1
    num(2) = find_name(handles.species_name1,'CO2');
    num(7) = find_name(handles.species_name2,'CO2');
else
    num(2) = -0.5;
end

if get(handles.OH,'value')==1
    num(3) = find_name(handles.species_name1,'OH');
    num(8) = find_name(handles.species_name2,'OH');
else
    num(3) = -0.5;    
end

if get(handles.CO,'value')==1
    num(4) = find_name(handles.species_name1,'CO');
    num(9) = find_name(handles.species_name2,'CO');
else
    num(4) = -0.5;    
end

if get(handles.H3O,'value')==1
    num(5) = find_name(handles.species_name1,'H3O+');
    num(10) = find_name(handles.species_name2,'H3O+');
else
    num(5) = -0.5;    
end


if get(handles.add_complex,'value')==1
    complex_label =1;
else
    complex_label =0;
end

% check whether all kinds of species in the network
for i =1:5
    if num(i) ==-1
        s = sprintf('No such kind of species in current chemistry network!');
        h = msgbox(s,'Warning Message','warn');
        return;
    end
end


k=1;
for i=1:5
    if num(i)>0
        axes(handles.axes1)
        semilogy(handles.z,handles.sigma1(:,num(i))./handles.abundance_H1(:), ...
            'color',handles.mymap(50*i,:),'LineWidth',2.5);
        legend_name(k) = handles.species_name1(num(i));
        k = k+1;
        hold on;
        if complex_label ==1
            semilogy(handles.z,handles.sigma4(:,num(i+5))./handles.abundance_H4(:), ...
                'color',handles.mymap(50*i,:),'LineStyle',':','LineWidth',2.5);
            legend_name(k) =strcat(handles.species_name2(num(i+5)),'-cplx');
            k = k +1;
        end
        
        axes(handles.axes2)
        semilogy(handles.z,handles.sigma2(:,num(i))./handles.abundance_H2(:), ...
            'color',handles.mymap(50*i,:),'LineWidth',2.5);
        hold on;
        if complex_label ==1
            semilogy(handles.z,handles.sigma5(:,num(i+5))./handles.abundance_H5(:), ...
                'color',handles.mymap(50*i,:),'LineStyle',':','LineWidth',2.5);
        end
        
        axes(handles.axes3)
        semilogy(handles.z,handles.sigma3(:,num(i))./handles.abundance_H3(:), ...
            'color',handles.mymap(50*i,:),'LineWidth',2.5);
        hold on;
        if complex_label ==1
            semilogy(handles.z,handles.sigma6(:,num(i+5))./handles.abundance_H6(:), ...
                'color',handles.mymap(50*i,:),'LineStyle',':','LineWidth',2.5);
        end
        
    end
end

axes(handles.axes1)
xlabel('Scale Height');
ylabel('Relative Abundance /[H]');
handles.legendObject1 = legend(legend_name,'Location','Southeast');
set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
hold off;
axes(handles.axes2)
xlabel('Scale Height');
ylabel('Relative Abundance /[H]');
handles.legendObject2= legend(legend_name,'Location','Southeast');
set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
hold off;
axes(handles.axes3)
xlabel('Scale Height');
ylabel('Relative Abundance /[H]');
handles.legendObject3 = legend(legend_name,'Location','Southeast');
set(gca,'linewidth',2,'fontsize',15,'fontname','Times');
hold off;
guidata(hObject,handles);


function vertical_range_Callback(hObject, eventdata, handles)
% hObject    handle to vertical_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vertical_range as text
%        str2double(get(hObject,'String')) returns contents of vertical_range as a double
str = get(hObject,'string');
str = regexp(str,'\s+','split');
if length(str) ~=2
    s = sprintf('Please input a right vertical range!');
    h = msgbox(s,'Warning Message','warn');
    return;
end
for i =1:2
    num(i) = str2double(str(i));
end
if num(1)>num(2)
    s = sprintf('Please set a right vertical range!');
    h = msgbox(s,'Warning Message','warn');
    return;
end
axes(handles.axes1)
axis([0,handles.height,num(1), num(2)]);
axes(handles.axes2)
axis([0,handles.height,num(1),num(2)]);
axes(handles.axes3)
axis([0,handles.height,num(1),num(2)]);

% --- Executes during object creation, after setting all properties.
function vertical_range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vertical_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CO.
function CO_Callback(hObject, eventdata, handles)
% hObject    handle to CO (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CO


% --- Executes on button press in add_complex.
function add_complex_Callback(hObject, eventdata, handles)
% hObject    handle to add_complex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add_complex


% --- Executes on button press in add_reduce_network.
function add_reduce_network_Callback(hObject, eventdata, handles)
% hObject    handle to add_reduce_network (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add_reduce_network



function NetLabel_Callback(hObject, eventdata, handles)
% hObject    handle to NetLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NetLabel as text
%        str2double(get(hObject,'String')) returns contents of NetLabel as a double

str = get(hObject,'string');
str = regexp(str,'\s+','split');
if length(str) ~=2
    s = sprintf('Please input right number of Chemical network');
    h = msgbox(s,'Warning Message','warn');
    return;
end
netname1 = str(1);
netname2 = str(2);
handles.netname1 = netname1;
handles.netname2 = netname2;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function NetLabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NetLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Vertical_num_Callback(hObject, eventdata, handles)
% hObject    handle to Vertical_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vertical_num as text
%        str2double(get(hObject,'String')) returns contents of Vertical_num as a double
str = get(hObject,'string');
str = regexp(str,'\s+','split');
if length(str) ~=1
    s = sprintf('Please input right vertical number');
    h = msgbox(s,'Warning Message','warn');
    return;
end
length1= str2double(str(1));
handles.length = length1;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Vertical_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vertical_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Scale_height_Callback(hObject, eventdata, handles)
% hObject    handle to Scale_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Scale_height as text
%        str2double(get(hObject,'String')) returns contents of Scale_height as a double
str = get(hObject,'string');
str = regexp(str,'\s+','split');
if length(str) ~=1
    s = sprintf('Please input right vertical height!');
    h = msgbox(s,'Warning Message','warn');
    return;
end
height= str2double(str(1));
handles.height= height;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Scale_height_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scale_height (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ReadData.
function ReadData_Callback(hObject, eventdata, handles)
% hObject    handle to ReadData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = handles.netname1;
[sigma1,species_name1] = readfile(str,handles.length,1);
handles.sigma1 =sigma1;
handles.species_name1 = species_name1;

[sigma2,species_name2] = readfile(str,handles.length,10);
handles.sigma2 =sigma2;

[sigma3,species_name3] = readfile(str,handles.length,100);
handles.sigma3 =sigma3;

str = handles.netname2;
[sigma4,species_name4] = readfile(str,handles.length,1);
handles.sigma4 =sigma4;
handles.species_name2 = species_name4;

[sigma5,species_name5] = readfile(str,handles.length,10);
handles.sigma5 =sigma5;

[sigma6,species_name6] = readfile(str,handles.length,100);
handles.sigma6 =sigma6;



% readfile(network, num_of_file, AU)
% 0 indicate simple network, 1 for complex network

%%% Z is the height /Unit in scale height
z=linspace(0,handles.height,handles.length);
handles.z = z;

% zetamax = 1e-5;
% zetamin = 1e-19;
% nzeta = 50;
% for i=0:handles.length-1
%     z(i+1) = zetamin*exp(i*log(zetamax/zetamin)/(handles.length-1));
% end
% 
% handles.z = log10(z);



% Get the total abundance of Hydrogen
for k=1:handles.length
    % hydrogen abundance in simple chemical network
    abundance_H1(k) = Find_element_abundance(species_name1,sigma1,k,'H');
    abundance_H2(k) = Find_element_abundance(species_name1,sigma2,k,'H');
    abundance_H3(k) = Find_element_abundance(species_name1,sigma3,k,'H');

    %hydrogen abundance in complex chemical network.
    abundance_H4(k) = Find_element_abundance(species_name4,sigma4,k,'H');
    abundance_H5(k) = Find_element_abundance(species_name4,sigma5,k,'H');
    abundance_H6(k) = Find_element_abundance(species_name4,sigma6,k,'H');
end
handles.abundance_H1 = abundance_H1;
handles.abundance_H2 = abundance_H2;
handles.abundance_H3 = abundance_H3;
handles.abundance_H4 = abundance_H4;
handles.abundance_H5 = abundance_H5;
handles.abundance_H6 = abundance_H6;

%%%%%%%% Set colorbar %%%%%%%%%%%%%%%%%%%%%%%%%
mymap = jet(256);
handles.mymap =mymap;

 guidata(hObject,handles);


% --- Executes on button press in SavePlot.
function SavePlot_Callback(hObject, eventdata, handles)
% hObject    handle to SavePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[filename, pathname] = uiputfile({ '*.png','Enhanced Meta File (*.png)';...
'*.bmp','Bitmap (*.bmp)'; '*.fig','Figure (*.fig)'}, ... 
'Save picture as','default');

% ??????????????
if isequal(filename,0) || isequal(pathname,0)
return
end
% ??????figure
newFig = figure;

% ??axes????????
axes_units = get(handles.axes1,'Units');
axes_pos = get(handles.axes1,'Position');

% ???????????figure?
axesObject2 = copyobj(handles.axes1,newFig);

% ?????????figure???????
set(axesObject2,'Units',axes_units);
set(axesObject2,'Position',[15 5 axes_pos(3) axes_pos(4)]);

legendObject = handles.legendObject;
% ??legendObject??????
%if (exist('legendObject'))

% ??legend??????
legend_units = get(legendObject,'Units');
legend_pos = get(legendObject,'Position');

% ?legend?????figure
legendObject2 = copyobj([legendObject],newFig);

% ????legend??
set(legendObject2,'Units',legend_units);
%set(legendObject2,'Position',[15 5 legend_pos(3) legend_pos(4)]);
set(legendObject2,'Position',[15-axes_pos(1)+legend_pos(1) 5-axes_pos(2)+legend_pos(2) legend_pos(3) legend_pos(4)] );
%end


% ???figure??????
set(newFig,'Units',axes_units);
set(newFig,'Position',[15 5 axes_pos(3)+30 axes_pos(4)+10]);

% ????
saveas(newFig,fullfile(pathname, filename)) 
close(newFig)

% --- Executes on selection change in SaveList.
function SaveList_Callback(hObject, eventdata, handles)
% hObject    handle to SaveList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SaveList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SaveList

value = get(hObject,'value');
handles.PlotName = value;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function SaveList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaveList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
