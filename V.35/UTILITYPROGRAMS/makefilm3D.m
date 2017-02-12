
% Resize the figure window to size of movie required.

% Record the size of the plot window:

% Within a loop, plot each picture and save to MATLAB movie matrix:
j=0;
aviobj = avifile('FILM.avi')
for j=0:48
	i = 48 - j;
	adderaruta = 0;
 	if ( i == 48 ) 
		adderaruta = 1;
		open ALANIN.48.fig
	end

 	if ( i == 46 )
	        adderaruta = 1;	
		open ALANIN.46.fig
	end
 	if ( i == 44 ) 
		adderaruta = 1;
		open ALANIN.44.fig
	end
 	if ( i == 42 ) 
		adderaruta = 1;
		open ALANIN.42.fig
	end
 	if ( i == 40 ) 
		adderaruta = 1;
		open ALANIN.40.fig
	end
 	if ( i == 38 ) 
		adderaruta = 1;
		open ALANIN.38.fig
	end
 	if ( i == 30 ) 
		adderaruta = 1;
		open ALANIN.30.fig
	end
 	if ( i == 28 ) 
		adderaruta = 1;
		open ALANIN.28.fig
	end
 	if ( i == 22 ) 
		adderaruta = 1;
		open ALANIN.22.fig
	end
 	if ( i == 20 ) 
		adderaruta = 1;
		open ALANIN.20.fig
	end
 	if ( i == 18 ) 
		adderaruta = 1;
		open ALANIN.18.fig
	end
 	if ( i == 16 ) 
		adderaruta = 1;
		open ALANIN.16.fig
	end
 	if ( i == 14 ) 
		adderaruta = 1;
		open ALANIN.14.fig
	end
 	if ( i == 12 ) 
		adderaruta = 1;
		open ALANIN.12.fig
	end
 	if ( i == 6 ) 
		adderaruta = 1;
		open ALANIN.6.fig
	end
 	if ( i == 4 ) 
		adderaruta = 1;
		open ALANIN.4.fig
	end
 	if ( i == 0 ) 
		adderaruta = 1;
		open ALANIN.0.fig
	end
	if ( adderaruta == 1 )
		axis([-6 8 -5 7 -6 4 ]);
		scrsz = get(0,'ScreenSize');
		%figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
		fig1=figure(1);
		winsize = get(fig1,'Position');
		set(fig1,'NextPlot','replacechildren')
		F = getframe(fig1);
		for k=1:16
			aviobj = addframe(aviobj,F);
		end
		close(fig1);
	end
end
aviobj = close(aviobj);
%close(fig1);
%This procedure creates a movie stored in a special format that is only readable in MATLAB. The first thing you will want to do is to play the movie:
