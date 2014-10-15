---
layout: lesson
root: ../..
title: HTML
---
<div class="objectives" markdown="1">

#### Learning Objectives
*   Explain the difference between text, elements, and tags.
*   Explain the difference between a model and a view, and correctly identify instances of each.
*   Write correctly-formatted HTML (using escape sequences for special characters).
*   Identify and fix improperly-nested HTML.
*   Explain what element attributes are, and what they are for.
*   Write HTML that uses attributes to alter a document's appearance.
*   Explain when to use attributes rather than nested elements.
*   Write correctly-formatted HTML pages containing lists, tables, images, and links.
*   Add correctly-formatted metadata to the head of an HTML page.

</div>

A basic HTML [document](../../gloss.html#document)
contains [text](../../gloss.html#text)
and [elements](../../gloss.html#element).
(The full specification allows for many other things
with names like "external entity references" and "processing instructions",
but we'll ignore them.)
The text in a document is just characters,
and as far as HTML is concerned,
it has no intrinsic meaning:
"Feynman" is just seven characters,
not a person.

Elements are [metadata](../../gloss.html#metadata)
that describe the meaning of the document's content.
For example,
one element might signal a heading,
while another might indicate that something is a cross-reference.

Elements are written using [tags](../../gloss.html#tag-xml),
which must be enclosed in angle brackets `<...>`.
For example, `<cite>` is used to mark the start of a citation,
and `</cite>` is used to mark its end.
Elements must be properly nested:
if an element called `inner` begins inside an element called `outer`,
`inner` must end before `outer` ends.
This means that `<outer>...<inner>...</inner></outer>` is legal HTML,
but `<outer>...<inner>...</outer></inner>` is not.

Here are some commonly-used HTML tags:

<table>
  <tr><th>Tag</th>               <th>Usage</th></tr>
  <tr><td>`html`</td> <td>Root element of entire HTML document.</td></tr>
  <tr><td>`body`</td> <td>Body of page (i.e., visible content).</td></tr>
  <tr><td>`h1`</td>   <td>Top-level heading.</td></tr>
  <tr><td>`h2, h3, ...`</td>   <td>Lower-level heading.</td></tr>
  <tr><td>`p`</td>    <td>Paragraph.</td></tr>
  <tr><td>`em`</td>   <td>Emphasized text.</td></tr>
</table>

Finally,
every well-formed document started with a `DOCTYPE` declaration,
which looks like:

<pre>
<!DOCTYPE html>
</pre>

This tells programs what kind of elements are allowed to appear in the document:
`html` (by far the most common case),
`math` for MathML,
and so on.
Here is a simple HTML document that uses everything we've seen so far:

~~~
<!DOCTYPE html>
<html>
<body>
<h1>Dimorphism</h1>
<p>Occurring or existing in two different <em>forms</em>.</p>
</body>
</html>
~~~

A web browser like Firefox might present this document like this:

FIXME: screenshot

Other devices will display it differently.
A phone,
for example,
might use a different background color for the heading,
while a screen reader for people with visual disabilities
would read the text aloud.

There are a couple of other formatting rules we need to know
in order to create and understand documents.
The first is that browsers and other programs will usually ignore whitespace in documents,
so this:

~~~
<!DOCTYPE html>
<html>
  <body>
    <h1>Dimorphism</h1>
    <p>Occurring or existing in two different <em>forms</em>.</p>
  </body>
</html>
~~~

will be displayed in exactly the same way as the first version.
Most people find that the indentation in the second makes it easier to read.

Second,
we must use [escape sequences](../../gloss.html#escape-sequence)
to represent the special characters `<` and `>`
for the same reason that we have to use `\&quot;`
inside a double-quoted string in a program.
In HTML and XML,
an escape sequence is an ampersand '&amp;'
followed by the abbreviated name of the character
(such as 'amp' for "ampersand")
and a semi-colon.
The four most common escape sequences are:

<table>
  <tr><th>Sequence</th> <th>Character</th></tr>
  <tr><td>&amp;lt;</td> <td><</td></tr>
  <tr><td>&amp;gt;</td> <td>></td></tr>
  <tr><td>&amp;quot;</td> <td>'</td></tr>
  <tr><td>&amp;amp;</td> <td>&</td></tr>
</table>

One final formatting rule is that
every document must have a single [root element](../../gloss.html#root-element),
i.e., a single element must enclose everything else.
A document like this is therefore not strictly legal:

~~~
<h1>Dimorphism</h1>
<p>Occurring or existing in two different <em>forms</em>.</p>
~~~

because it has two top-level elements
(the `h1` and the `p`).
Most browsers will render it correctly,
though,
since they're designed to accommodate improperly-formatted HTML,
but most other programs will complain.

> #### Beautiful Soup
>
> There are a lot of incorrectly-formatted HTML pages out there.
> To deal with them,
> people have written libraries like <a href="http://www.crummy.com/software/BeautifulSoup/">Beautiful Soup</a>,
> which does its best to turn real-world HTML into something that
> a run-of-the-mill program can handle.
> It almost always gets things right,
> but sticking to the standard makes life a lot easier for everyone.

Elements can be customized by giving them [attributes](../../gloss.html#attribute).
These are name/value pairs enclosed in the opening tag like this:

~~~
<h1 align="center">A Centered Heading</h1>
~~~

or:

~~~
<p class="disclaimer">This planet provided as-is.</p>
~~~

Any particular attribute name may appear at most once in any element,
so `<p align="left" align="right">...</p>` is illegal.
Attributes' values <em>must</em> be in quotes in XML and older dialects of HTML;
HTML5 allows single-word values to be unquoted,
but quoting is still recommended.

Another important feature of attributes is that they are unordered.
They have to be <em>written</em> in some order,
but as far as the rules of HTML are concerned:

~~~
<p align="center" class="disclaimer">This web page is made from 100% recycled pixels.</p>
~~~

and:

~~~
<p class="disclaimer" align="center">This web page is made from 100% recycled pixels.</p>
~~~

mean the same thing.

> #### HTML and Version Control
> 
> FIXME

When should we use attributes, and when should we nest elements?
As a general rule,
we should use attributes when:

*   each value can occur at most once for any element;
*   the order of the values doesn't matter; and
*   those values have no internal structure, i.e., we will never need to parse an attribute's value in order to understand it.

In all other cases, we should use nested elements.
However, many widely-used XML formats break these rules
in order to make it easier for people to write XML by hand.
For example,
in the Scalable Vector Graphics (SVG) format used to describe images as XML,
we would define a rectangle as follows:

~~~
<rect width="300" height="100" style="fill:rgb(0,0,255); stroke-width:1; stroke:rgb(0,0,0)"/>
~~~

In order to understand the `style` attribute,
a program has to somehow know to split it on semicolons,
and then to split each piece on colons.
This means that a generic program for reading XML
can't extract all the information that's in SVG,
which partly defeats the purpose of using XML in the first place.

Web pages can contain a lot more than just headings and italics.
To start with,
HTML provides two kinds of lists:
`ul` to mark an unordered (bulleted) list,
and `ol` for an ordered (numbered) one.
Items inside either kind of list must be wrapped in `li` elements:

~~~
<!DOCTYPE html>
<html>
<body>
  <ul>
    <li>A. Binet
      <ol>
        <li>H. Ebbinghaus</li>
        <li>W. Wundt</li>
      </ol>
    </li>
    <li>C. S. Pierce
      <ol>
        <li>W. Wundt</li>
      </ol>
    </li>
</body>
</html>
~~~

FIXME: image

Note how elements are nested:
since the ordered lists "belong" to the unordered list items above them,
they are inside those items' `<li>...</li>` tags.

HTML also provides tables, but they are awkward to use:
tables are naturally two-dimensional,
but text is one-dimensional.
This is exactly like the problem of representing a two-dimensional array in memory,
and we solve it in the same way:
by writing down the rows,
and the columns within each row,
in a fixed order.
The `table` element marks the table itself;
within that,
each row is wrapped in `tr` (for "table row"),
and within those,
column items are wrapped in `th` (for "table heading")
or `td` (for "table data"):

~~~
<!DOCTYPE html>
<html>
<body>
  <table>
    <tr>
      <th></th>
      <th>A. Binet</th>
      <th>C. S. Pierce</th>
    </tr>
    <tr>
      <th>H. Ebbinghaus</th>
      <td>88%</td>
      <td>NA</td>
    </tr>
    <tr>
      <th>W. Wundt</th>
      <td>29%</td>
      <td>45%</td>
    </tr>
  </table>
</body>
</html>
~~~

FIXME: table

> #### Tables, Layout, and CSS
> 
> Tables are sometimes used to do multi-column layout,
> as well as for tabular data,
> but this is a bad idea.
> To understand why,
> consider two other HTML tags:
> `i`, meaning "italics",
> and `em`, meaning "emphasis".
> The former directly controls how text is displayed,
> but by doing so,
> it breaks the separation between model and view that is the heart of markup's usefulness.
> Without understanding the text that has been italicized,
> a program cannot understand whether it is meant to indicate someone shouting,
> the definition of a new term,
> or the title of a book.
> The `em` tag, on the other hand, has exactly one meaning,
> and that meaning is different from the meaning of `dfn` (a definition)
> or `cite` (a citation).
> 
> Conscientious authors use [Cascading Style Sheets](../../gloss.html#css) (or CSS)
> to describe how they want pages to appear,
> and only use `table` elements for actual tables.
> CSS is beyond the scope of this lesson,
> but is described briefly in <a href="extras.html#s:web:Css">The Appendix</A>.

HTML pages can also contain images.
(In fact,
the World Wide Web didn't really take off until
the Mosaic browser allowed people to mix images with text.)
The word "contain" is misleading, though:
HTML documents can only contain text,
so we cannot store an image "in" a page.
Instead,
we must put it in some other file,
and insert a reference to that file in the HTML using the `img` tag.
Its `src` attribute specifies where to find the image file;
this can be a path to a file on the same host as the web page,
or a URL for something stored elsewhere.
For example,
when a browser displays this:

~~~
<!DOCTYPE html>
<html>
  <body>
    <p>My daughter's first online chat:</p>
    <img src="img/maddie-chat.jpg"/>
    <p>but probably not her last.</p>
  </body>
</html>
~~~

it looks for the file `maddie-chat.jpg`
in the same directory as the HTML file:

FIXME: screenshot

Notice,
by the way,
that the `img` element is written as
`<img.../>`,
i.e.,
with a trailing slash inside the `<>`
rather than with a separate closing tag.
This makes sense because the element doesn't contain any text:
the content is referred to by its `src` attribute.
Any element that doesn't contain anything
can be written using this short form.

Notice also that images don't have to be in the same directory as the pages that refer to them:
the browser will follow relative paths (ones that don't start with `/`).

> #### It's Always Interpreted
> 
> FIXME: Absolute paths - the path is <em>always</em> interpreted (web browser config).

Whenever we refer to an image,
we should use the `img` tag's `alt` attribute
to provide a title or description of the image.
This is what screen readers for people with visual handicaps will say aloud to "display" the image;
it's also what search engines rely on,
since they can't "see" the image either.
Adding this to our previous example gives:

~~~
<!DOCTYPE html>
<html>
<body>
  <p>My daughter's first online chat:</p>
  <img src="maddie-chat.jpg" alt="Madeleine's first online chat"/>
  <p>but probably not her last.</p>
</body>
</html>
~~~

The other element --- the one that makes HTML pages "hypertext" --- is `a`.
Whatever is inside the element is displayed and highlighted for clicking;
this is usually a few words of text,
but it can be an entire paragraph or an image.

The `a` element's `href` attribute
specifies what the link is pointing at;
as with images,
this can be either a local filename or a URL.
For example,
we can create a listing of the examples we've written so far like this:

~~~
<!DOCTYPE html>
<html>
<body>
  <p>
    Simple HTML examples for
    <a href="http://software-carpentry.org">Software Carpentry</a>.
  </p>
  <ol>
    <li><a href="very-simple.html">a very simple page</a></li>
    <li><a href="hide-paragraph.html">hiding paragraphs</a></li>
    <li><a href="nested-lists.html">nested lists</a></li>
    <li><a href="simple-table.html">a simple table</a></li>
    <li><a href="simple-image.html">a simple image</a></li>
  </ol>
</body>
</html>
~~~

FIXME: image

The hyperlink element is called `a` because
it can also used to create [anchors](../../gloss.html#anchor) in documents
by giving them a `name` attribute instead of an `href`.
An anchor is simply a location in a document that can be linked to.
For example,
suppose we formatted the Feynman quotation given earlier like this:

~~~
<blockquote>
As a by-product of this same view, I received a telephone call one day
at the graduate college at <a name="princeton-university">Princeton</a>
from Professor Wheeler, in which he said,
"Feynman, I know why all electrons have the same charge and the same mass."
"Why?"
"Because, they are all the same electron!"
</blockquote>
~~~

If this quotation was in a file called `quote.html`,
we could then create a hyperlink directly to the mention of Princeton
using `<a&nbsp;href="quote.html#princeton-university">`.
The `#` in the `href`'s value separates the path to the document
from the anchor we're linking to.
Inside `quote.html` itself,
we could link to that same location simply using
`<a&nbsp;href="#pu">`.

Using the `a` element for both links and targets was poor design:
programs are simpler to write if each element has one purpose, and one alone.
A better way to create anchors is to add an `id` attribute
to some other element.
For example,
if we wanted to be able to link to the quotation itself,
we could write:

~~~
<blockquote id="wheeler-electron-quote">
As a by-product of this same view, I received a telephone call one day
at the graduate college at <span id="princeton-university">Princeton</span>
from Professor Wheeler, in which he said,
"Feynman, I know why all electrons have the same charge and the same mass."
"Why?"
"Because, they are all the same electron!"
</blockquote>
~~~

and then refer to `quote.html#wheeler-electron-quote`
or `quote.html#princeton-university`.

Finally,
well-written HTML pages have a `head` element as well as a `body`.
The head isn't displayed;
instead,
it's used to store metadata about the page as a whole.
The most common element inside `head` is `title`,
which,
as its name suggests,
gives the page's title.
(This is usually displayed in the browser's title bar.)
Another common item in the head is `meta`,
whose two attributes `name` and `content`
let authors add arbitrary information to their pages.
If we add these to the web page we wrote earlier,
we might have:

~~~
<!DOCTYPE html>
<html>
  <head>
    <title>Dimorphism Defined<title>
    <meta name="author" content="Alan Turing"/>
    <meta name="institution" content="Euphoric State University"/>
  </head>
  <body>
    <h1>Dimorphism</h1>
    <p>Occurring or existing in two different <em>forms</em>.</p>
  </body>
</html>
~~~

Well-written pages also use comments (just like code),
which start with `<!--` and end with `-->`.

> #### Hiding Content
> 
> Commenting out part of a page does <em>not</em> hide the content
> from people who really want to see it:
> while a browser won't display what's inside a comment,
> it's still in the page,
> and anyone who uses "View Source" can read it.
> For example,
> if you are looking at this page in a web browser right now,
> try viewing the source
> and searching for the word "Surprise".
> 
> <!-- Surprise: this isn't displayed by the browser, but is still in the document. -->
> 
> If you really don't want people to be able to read something,
> the only safe thing to do is to keep it off the web.

<div class="keypoints" markdown="1">

#### Summary
*   HTML documents contain elements and text.
*   Elements are represented using tags.
*   Different devices may display HTML differently.
*   Every document must have a single root element.
*   Special characters must be written using escape sequences beginning with &.
*   Elements can be customized by adding key-value pairs called attributes.
*   An element's attributes must be unique, and are unordered.
*   Attribute values should not have any internal structure.
*   Put metadata in `meta` elements in a page's `head` element.
*   Use `ul` for unordered lists and `ol` for ordered lists.
*   Add comments to pages using `<!--` and `-->`.
*   Use `table` for tables, with `tr` for rows and `td` for values.
*   Use `img` for images.
*   Use `a` to create hyperlinks.
*   Give elements a unique `id` attribute to link to it.

</div>
