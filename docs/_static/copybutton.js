/*
 * @name    Copybutton
 * @desc    Add a button to the code examples in the doc to show/hide the prompts and output.
 * @author  Ezio Melotti <ezio.melotti@gmail.com>
 * @author  Georg Brandl <georg@python.org>
 * @date    Sun, 04 Dec 2011 11:51:21 +0100 (2011-12-04)
 * @license Python Software Foundation License (PSFL) <https://docs.python.org/2/license.html>
 * @url     <http://hg.python.org/cpython/file/0436ef8be253/Doc/tools/sphinxext/static/copybutton.js>
 * @home    <http://hg.python.org/cpython>
 * 
 * Version history:
 *   (2011-12-04) Georg Brandl <http://hg.python.org/cpython/rev/0436ef8be253>
 *   (2011-10-30) Ezio Melotti <http://hg.python.org/cpython/rev/18bbfed9aafa>
 * 
 */


$(document).ready(function() {
    /* Add a [>>>] button on the top-right corner of code samples to hide
     * the >>> and ... prompts and the output and thus make the code
     * copyable. */
    var div = $('.highlight-python .highlight,' +
                '.highlight-default .highlight,' +
                '.highlight-python3 .highlight')
    var pre = div.find('pre');

    // get the styles from the current theme
    pre.parent().parent().css('position', 'relative');
    var hide_text = 'Hide the prompts and output';
    var show_text = 'Show the prompts and output';
    var border_width = pre.css('border-top-width');
    var border_style = pre.css('border-top-style');
    var border_color = pre.css('border-top-color');
    var button_styles = {
        'cursor':'pointer', 'position': 'absolute', 'top': '0', 'right': '0',
        'border-color': border_color, 'border-style': border_style,
        'border-width': border_width, 'color': border_color, 'text-size': '75%',
        'font-family': 'monospace', 'padding-left': '0.2em', 'padding-right': '0.2em'
    }

    // create and add the button to all the code blocks that contain >>>
    div.each(function(index) {
        var jthis = $(this);
        if (jthis.find('.gp').length > 0) {
            var button = $('<span class="copybutton">&gt;&gt;&gt;</span>');
            button.css(button_styles)
            button.attr('title', hide_text);
            jthis.prepend(button);
        }
        // tracebacks (.gt) contain bare text elements that need to be
        // wrapped in a span to work with .nextUntil() (see later)
        jthis.find('pre:has(.gt)').contents().filter(function() {
            return ((this.nodeType == 3) && (this.data.trim().length > 0));
        }).wrap('<span>');
    });

    // define the behavior of the button when it's clicked
    $('.copybutton').click(function() {
        var button = $(this);
        if (!button.hasClass('copybutton-js-clicked')) {
            button.addClass('copybutton-js-clicked');
            button.parent().find('.go, .gp, .gt').hide();
            button.next('pre').find('.gt').nextUntil('.gp, .go').css('visibility', 'hidden');
            button.css('text-decoration', 'line-through');
            button.attr('title', show_text);
        } else {
            button.removeClass('copybutton-js-clicked');
            button.parent().find('.go, .gp, .gt').show();
            button.next('pre').find('.gt').nextUntil('.gp, .go').css('visibility', 'visible');
            button.css('text-decoration', 'none');
            button.attr('title', hide_text);
        }
    });
});

