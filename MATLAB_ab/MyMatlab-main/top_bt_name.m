function top_bt_name(titl,txtup,txtbot,postit,fsz)

% function top_bt_name(titl,txtup,txtbot,postit,fsz)
% creates subscripts on the plot
% titl - title of the figure
% txtup - text string that will be placed in the upper right coner
% txtbot - test string that will be placed in the lower left coner
% postit - the coordinates for the title, if they are not specified-> [], 
% fsz - fontsize of the title;
% program uses default xtit, ytit
% date is placed in the left lower coner

	if (isempty(postit))
	  xtit=0.5;
	  ytit=1.05;
	else
	  xtit=postit(1);
	  ytit=postit(2);
	end

  if (isempty(fsz)), fsz=12; end;

% do bottom program name
	h=axes('Units', 'normalized','Position',[.1 .1 .8 .8]);
	axes(h);
	axis off
	h=gca;
 	set(h,'Clipping','off')
	set(h,'FontUnits','points')
	h=text(0., -0.1, [datestr(now,2)],'HorizontalALignment','left',...
        		'FontSize',6,'Units','normalized');

	h=text(xtit,ytit,titl,...
		'HorizontalAlignment','center',...
		'FontSize',fsz,...
		'Units','normalized'); 

  	h=text(1.1,-0.1,txtbot,...
      		'HorizontalAlignment','right',...
      		'FontSize',6,...
      		'Units','normalized');

 	h=text(1.,1.10,txtup,...
      		'HorizontalAlignment','right',...
      		'FontSize',6,...
      		'Units','normalized');
 


