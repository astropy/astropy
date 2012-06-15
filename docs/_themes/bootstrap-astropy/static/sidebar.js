/*
 * sidebar.js
 * ~~~~~~~~~~
 *
 * This script makes the Sphinx sidebar collapsible.
 *
 * .sphinxsidebar contains .sphinxsidebarwrapper.  This script adds
 * in .sphixsidebar, after .sphinxsidebarwrapper, the #sidebarbutton
 * used to collapse and expand the sidebar.
 *
 * When the sidebar is collapsed the .sphinxsidebarwrapper is hidden
 * and the width of the sidebar and the margin-left of the document
 * are decreased. When the sidebar is expanded the opposite happens.
 * This script saves a per-browser/per-session cookie used to
 * remember the position of the sidebar among the pages.
 * Once the browser is closed the cookie is deleted and the position
 * reset to the default (expanded).
 *
 * :copyright: Copyright 2007-2011 by the Sphinx team, see AUTHORS.
 * :license: BSD, see LICENSE for details.
 *
 */

$(function() {
  // global elements used by the functions.
  // the 'sidebarbutton' element is defined as global after its
  // creation, in the add_sidebar_button function
  var bodywrapper = $('.bodywrapper');
  var sidebar = $('.sphinxsidebar');
  var sidebarwrapper = $('.sphinxsidebarwrapper');

  // for some reason, the document has no sidebar; do not run into errors
  if (!sidebar.length) return;

  // original margin-left of the bodywrapper and width of the sidebar
  // with the sidebar expanded
  var bw_margin_expanded = bodywrapper.css('margin-left');
  var ssb_width_expanded = sidebar.width();

  // margin-left of the bodywrapper and width of the sidebar
  // with the sidebar collapsed
  var bw_margin_collapsed = 12;
  var ssb_width_collapsed = 12;

  // custom colors
  var dark_color = '#404040';
  var light_color = '#505050';

  function sidebar_is_collapsed() {
    return sidebarwrapper.is(':not(:visible)');
  }

  function toggle_sidebar() {
    if (sidebar_is_collapsed())
      expand_sidebar();
    else
      collapse_sidebar();
  }

  function collapse_sidebar() {
    sidebarwrapper.hide();
    sidebar.css('width', ssb_width_collapsed);
    bodywrapper.css('margin-left', bw_margin_collapsed);
    sidebarbutton.css({
        'margin-left': '-1px',
        'height': bodywrapper.height(),
	'border-radius': '3px'
    });
    sidebarbutton.find('span').text('»');
    sidebarbutton.attr('title', _('Expand sidebar'));
    document.cookie = 'sidebar=collapsed';
  }

  function expand_sidebar() {
    bodywrapper.css('margin-left', bw_margin_expanded);
    sidebar.css('width', ssb_width_expanded);
    sidebarwrapper.show();
    sidebarbutton.css({
        'margin-left': ssb_width_expanded - 12,
        'height': bodywrapper.height(),
	'border-radius': '0px 3px 3px 0px'
    });
    sidebarbutton.find('span').text('«');
    sidebarbutton.attr('title', _('Collapse sidebar'));
    document.cookie = 'sidebar=expanded';
  }

  function add_sidebar_button() {
    sidebarwrapper.css({
        'float': 'left',
        'margin-right': '0',
        'width': ssb_width_expanded - 18
    });
    // create the button
    sidebar.append('<div id="sidebarbutton"><span>&laquo;</span></div>');
    var sidebarbutton = $('#sidebarbutton');

    // find the height of the viewport to center the '<<' in the page
    var viewport_height;
    if (window.innerHeight)
 	  viewport_height = window.innerHeight;
    else
	  viewport_height = $(window).height();
    var sidebar_offset = sidebar.offset().top;
    var sidebar_height = Math.max(bodywrapper.height(), sidebar.height());
    sidebarbutton.find('span').css({
        'font-family': '"Lucida Grande",Arial,sans-serif', 
	'display': 'block',
	'top': Math.min(viewport_height/2, sidebar_height/2 + sidebar_offset) - 10,
	'width': 12,
	'position': 'fixed',
	'text-align': 'center'
    });

    sidebarbutton.click(toggle_sidebar);
    sidebarbutton.attr('title', _('Collapse sidebar'));
    sidebarbutton.css({
        'color': '#FFFFFF',
	'background-color': light_color,
	'border': '1px solid ' + light_color,
        'border-radius': '0px 3px 3px 0px',
        'font-size': '1.2em',
        'cursor': 'pointer',
        'height': sidebar_height,
        'padding-top': '1px',
	'margin': '-1px',
        'margin-left': ssb_width_expanded - 12
    });

    sidebarbutton.hover(
      function () {
          $(this).css('background-color', dark_color);
      },
      function () {
          $(this).css('background-color', light_color);
      }
    );
  }

  function set_position_from_cookie() {
    if (!document.cookie)
      return;
    var items = document.cookie.split(';');
    for(var k=0; k<items.length; k++) {
      var key_val = items[k].split('=');
      var key = key_val[0];
      if (key == 'sidebar') {
        var value = key_val[1];
        if ((value == 'collapsed') && (!sidebar_is_collapsed()))
          collapse_sidebar();
        else if ((value == 'expanded') && (sidebar_is_collapsed()))
          expand_sidebar();
      }
    }
  }

  add_sidebar_button();
  var sidebarbutton = $('#sidebarbutton');
  set_position_from_cookie();
});
