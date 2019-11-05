%CROP_MARGINS Remove margins from axes object
%
%   CROP_MARGINS(ax) changes the Position property of the axes object ax
%   so that there are no margins around the axis except those for title and
%   labels.
%
%   CROP_MARGINS() operates on the current axes from gca().

function crop_margins(ax)
  if nargin == 0
    ax = gca();
  end
  ax.Units = 'normalized';

  % Get margins around axis required for title and labels.
  % Add 0.5% of the axes width or height to each margin.
  extra = 0.005;
  room_left   = extra + ax.TightInset(1);
  room_bottom = extra + ax.TightInset(2);
  room_right  = extra + ax.TightInset(3);
  room_top    = extra + ax.TightInset(4);

  % Change position of inner axes so that with title and labels it occupies the
  % complete outer axes.
  ax.Position = ...
  [room_left; room_bottom; 1-(room_left+room_right); 1-(room_bottom+room_top)];
end
