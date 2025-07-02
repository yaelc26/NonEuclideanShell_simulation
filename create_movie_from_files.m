function M=create_movie_from_files(v_output,f_input,f_output,hide)

n_nodes=926; n_faces=1765; start_at=0; %Don't Start at beginning to avoid mem issues

if nargin<3, nofo=true; cmax=1; else nofo=false; end;

t=dlmread(f_input,'\t',1,0); t=t(:,2:4)';
p=dlmread(v_output,'\t',n_nodes*start_at,0); frames=unique(p(:,5))'; alims=minmax(p(:,2:4)')'; alims=alims(:)'+eps*[-1 1 -1 1 -1 1];
if ~nofo, c=dlmread(f_output,'\t',n_faces*start_at,0); nc=size(c,2)-2; cmax=max(c(:,2:end-1)); end;

for f=frames
    m=p(:,5)==f;
    x=p(m,2); y=p(m,3); z=p(m,4);
    M(f-start_at).x=x(t+1); M(f-start_at).y=y(t+1); M(f-start_at).z=z(t+1);
    if ~nofo
        for cindex=1:nc
            col=c(c(:,nc+2)==f,cindex+1)';
            M(f-start_at).c{cindex}=col;
        end
    end
end
%%%%%%%
%M = M(end);
%%%%%%%
if nargin<4||~hide, myGUI(M,alims,nofo,cmax); end;

end

function myGUI(M,alims,nofo,cmax)
% Build GUI
figH = figure(                ...
    'Units'     , 'Pixels'    , ...
    'Toolbar'     , 'Figure'    );
axesH = axes(                 ...
    'Units'     , 'Pixels'    );
sliderH = uicontrol(          ...
    'Style'     , 'Slider'    , ...
    'Units'     , 'Pixels'    , ...
    'Min'       , eps           , ...
    'Max'       , length(M)     , ...
    'Value'     , eps        , ...
    'SliderStep'     , [1/length(M) 1/7]      , ...
    'Callback'  , @mySliderFcn);
buttonH = uicontrol(          ...
    'Style'     , 'Pushbutton', ...
    'Units'     , 'Pixels'    , ...
    'String'    , 'view Eb'     , ...
    'Callback'  , @myButtonFcn);
if nofo, set(buttonH,'Enable','off'); end;
% Plot
toggle=1;
plotSurface(1);
cbarH=colorbar;
cP=get(cbarH  , 'Position');
set(figH,'ResizeFcn' , @myResizeFcn);
myResizeFcn;

    function myResizeFcn(varargin)
        % Figure resize callback
        %   Adjust object positions so that they maintain appropriate
        %   proportions
        fP = get(figH, 'Position');
        set(sliderH, 'Position', [10       , 10, fP(3)-130, 25       ]);
        set(buttonH, 'Position', [fP(3)-110, 10, 100      , 25       ]);
        set(axesH  , 'Position', [50       , 85, fP(3)-100, fP(4)-130]);
        set(cbarH  , 'Position', cP);
    end
    function mySliderFcn(varargin)
        % Slider callback
        %   Modifies the coefficient
        fr = ceil(get(sliderH, 'Value'));
        plotSurface(fr);
    end
    function myButtonFcn(varargin)
        % Button callback
        fr=ceil(get(sliderH,'Value'));
        toggle=mod(toggle,2)+1;
        switch toggle
            case 1
                set(buttonH,'String','view Eb')
            case 2
                set(buttonH,'String','view Es')
        end
        plotSurface(fr);
    end
    function plotSurface(fr)
        % Plot equation with given coefficient
        axes(axesH); cla;
        if nofo
            patch(M(fr).x,M(fr).y,M(fr).z,M(fr).z,'EdgeColor','k','FaceColor',[0 0.5 0]);
            caxis([0 cmax(toggle)]); %material dull; lighting phong; lightangle(30,30); lightangle(-30,30); lightangle(0,-45);
        else
            patch(M(fr).x,M(fr).y,M(fr).z,M(fr).c{toggle},'EdgeColor','none');
            caxis([0 cmax(toggle)]); material dull;  lighting phong; lightangle(90+30,30); lightangle(90-30,30); lightangle(0,-45);
        end
        axis equal; axis(alims); axis off; %view(80,20);
        switch toggle
            case 1
                title(['Stretching Energy Density - frame ' num2str(fr) ' of ' num2str(length(M))]);
            case 2
                title(['Bending Energy Density - frame ' num2str(fr) ' of ' num2str(length(M))]);
        end
    end
end
