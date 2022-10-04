\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{TUD-report2020}[2013/06/22 v1.0 TU Delft thesis class]

\newif\if@dutch
\newif\if@english
\newif\if@nativefonts
%\newif\if@whitetitle
% \newif\if@whitelogo
\DeclareOption{dutch}{\@dutchtrue\@englishfalse}
\DeclareOption{english}{\@dutchfalse\@englishtrue}
\DeclareOption{nativefonts}{\@nativefontstrue}
%\DeclareOption{whitetitle}{\@whitetitletrue}
% \DeclareOption{whitelogo}{\@whitelogotrue}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ExecuteOptions{english}
\ProcessOptions\relax

\LoadClass[11pt]{book}

\RequirePackage{amsmath}
\RequirePackage{amssymb}
%% English is the default language, but this can be changed to Dutch by
%% specifying the 'dutch' option to the document class.
\if@dutch
    \RequirePackage[dutch]{babel}
\else
    %% Even for English reports, some sections might be in Dutch, for example on
    %% the title page.
    \RequirePackage[dutch,english]{babel}
\fi
\RequirePackage[nooneline,footnotesize]{caption}
\RequirePackage{fancyhdr}
\RequirePackage[flushmargin,hang]{footmisc}
\RequirePackage{ifpdf}
\RequirePackage{ifxetex}
\ifxetex
    \RequirePackage[xetex]{geometry}
    \RequirePackage[xetex]{graphicx}
    \RequirePackage[xetex,breaklinks]{hyperref}
    \RequirePackage{fontspec}
    \RequirePackage{xltxtra}
    \defaultfontfeatures{Ligatures=TeX}
\else
    \ifpdf
        \RequirePackage[pdftex]{geometry}
        \RequirePackage[pdftex]{graphicx}
        \RequirePackage[pdftex,breaklinks]{hyperref}
    \else
        \RequirePackage[dvips]{geometry}
        \RequirePackage[dvips]{graphicx}
        \RequirePackage[hypertex,breaklinks]{hyperref}
    \fi
\fi
\RequirePackage{metalogo}
\RequirePackage{tikz}
\RequirePackage[noindentafter]{titlesec}
\RequirePackage{titletoc}
\RequirePackage{xcolor}
\RequirePackage{xkeyval}
\RequirePackage[export]{adjustbox}

\@nativefontstrue
\newcommand*\tudrmfamily{\rmfamily}
\newcommand*\tudsffamily{\sffamily}


%% The native font is Utopia Regular with Fourier math. We use Latin Modern
%% as a fallback, since Utopia lacks support for sans-serif and monospace.
\RequirePackage[T1]{fontenc}
\RequirePackage[utf8]{inputenc}

\RequirePackage{csquotes}

% load fonts
% \RequirePackage{fourier}
\RequirePackage{sourcecodepro}
\RequirePackage{FiraSans}

% Load font for the text of the entire document
\RequirePackage{mathpazo}

\def\titlefont{\firabook}
\def\sectitlestyle{\rmfamily\fontshape{n}\bfseries}
\def\chaptitlestyle{\rmfamily\fontshape{n}\bfseries}
\def\titlestyle{\titlefont\fontshape{n}}
\def\headerstyle{\titlefont}

% Use roman (Georgia) in quote environment
\renewenvironment{quote}
               {\list{}{\rightmargin\leftmargin}\tudrmfamily\selectfont%
                \item\relax}
               {\endlist}

%% Define the TU Delft house style colors.
\definecolor{tudelft-cyan}{cmyk}{1,0,0,0}
\definecolor{tudelft-black}{cmyk}{0,0,0,1}
\definecolor{tudelft-white}{cmyk}{0,0,0,0}
\definecolor{tudelft-sea-green}{cmyk}{0.54,0,0.32,0}
\definecolor{tudelft-green}{cmyk}{1,0.15,0.4,0}
\definecolor{tudelft-dark-blue}{cmyk}{1,0.66,0,0.4}
\definecolor{tudelft-purple}{cmyk}{0.98,1,0,0.35}
\definecolor{tudelft-turquoise}{cmyk}{0.82,0,0.21,0.08}
\definecolor{tudelft-sky-blue}{cmyk}{0.45,0,0.06,0.06}
\definecolor{tudelft-lavendel}{cmyk}{0.45,0.2,0,0.07}
\definecolor{tudelft-orange}{cmyk}{0.02,0.56,0.84,0}
\definecolor{tudelft-warm-purple}{cmyk}{0.58,1,0,0.02}
\definecolor{tudelft-fuchsia}{cmyk}{0.19,1,0,0.19}
\definecolor{tudelft-bright-green}{cmyk}{0.36,0,1,0}
\definecolor{tudelft-yellow}{cmyk}{0.02,0,0.54,0}

\def\@pagecolor{tudelft-cyan}
\colorlet{title}{tudelft-cyan}

%% We use A4 paper with slightly smaller margins than the default (scale = 0.7).
\geometry{a4paper,hscale=0.75,vscale=0.8}

%% Redefine the title command to accept an optional subtitle.
\renewcommand*\title[2][tudelft-white]{%
    \def\@titlecolor{#1}%
    \def\@title{#2}%
    %% Add the title to the PDF meta data.
    \hypersetup{pdftitle=#2}%
}
\newcommand*\subtitle[2][tudelft-white]{%
%    \def\@subtitle{#1}%
    \def\@subtitlecolor{#1}%
    \def\@subtitle{#2}%
    %% Add the title to the PDF meta data.
%    \hypersetup{pdftitle=#3}%
}
%% Redefine the author command to add the name to the PDF meta data.
\renewcommand*\author[2][tudelft-white]{%
    \def\@authorcolor{#1}%
    \def\@author{#2}%
    \hypersetup{pdfauthor=#2}%
}
%% The coverimage command is used to specify the filename of the optional cover
%% image.
\newcommand*\coverimage[1]{\def\@cover@image{#1}}

\newcommand\covertext[2][tudelft-white]{%
    \def\@covertextcolor{#1}%
    \def\@cover@text{#2}}


%% The cover image is scaled to fill the width of the front cover.
\newlength\@cover@fronttextwidth
\newlength\@cover@fronttextheight
\newlength\@cover@splitboxtextwidth
\newlength\@cover@splittextheight

%% The margin around text boxes.
\newlength\@cover@margin
\setlength\@cover@margin{20pt}


\newcommand*\makecover[1][]{
    \setkeys{cover}{#1}
    \clearpage
    %% Create a cyan-colored empty page without margins.
    \newgeometry{margin=0pt}
    \pagecolor{tudelft-cyan}
    \thispagestyle{empty}
    %% We need the this to perform coordinate calculations in TikZ.
    \usetikzlibrary{calc}
    \begin{tikzpicture}[remember picture,overlay]
        %% If a cover image was specified, attach it to the top right of the
        %% front cover.
        \ifx\@cover@image\undefined\else
            \node at (current page.south)[anchor=south,inner sep=0pt]{
            \includegraphics[width=\textwidth]{\@cover@image}};
        \fi
        \node at (current page.south west)[anchor=south west,inner xsep=\@cover@margin,inner ysep=10pt]{
          \includegraphics{cover/logo_white}
        };
        %% Calculate the coordinate of the top left corner of the front cover.
        \coordinate (top left) at (current page.north west);
        \coordinate (bottom left) at (current page.south west);
        %% Calculate the coordinate of the corner where the front and back boxes
        %% meet.
        \setlength\@cover@fronttextwidth{\textwidth-4\@cover@margin}
        \setlength\@cover@fronttextheight{0.45\paperheight-\@cover@margin}
        \coordinate (tcorner) at ($(top left)+(\@cover@margin,-\@cover@margin)$);
        \node at (tcorner)[anchor=north west,inner sep=\@cover@margin]{
            \begin{minipage}[t][\@cover@fronttextheight]{\@cover@fronttextwidth}
                %% Print the title and optional subtitle at the top in white.
                {\titlefont\color{\@titlecolor}\fontsize{50}{50}\selectfont\@title}
                \ifx\@subtitle\undefined\else
                    \vfill
                    {\titlefont\color{\@subtitlecolor}\fontsize{52}{52}\selectfont\@subtitle}
                \fi
                %% Print the author.
                \vfill
                \titlefont\color{\@authorcolor}\fontsize{50}{50}\selectfont\@author
                \vfill
                \titlefont\color{\@covertextcolor}\fontsize{18}{18}\selectfont\@cover@text
            \end{minipage}
        };
    \end{tikzpicture}
    %% Restore the margins and turn the page white again.
    \restoregeometry
    \pagecolor{white}
}

%% Remove the header and page number on empty pages.
\def\cleardoublepage{%
    \clearpage%
    \if@twoside%
        \ifodd\c@page%
        \else%
            \thispagestyle{empty}%
            \vspace*{\fill}%
            \newpage%
        \fi%
    \fi%
}

%% Page style for title pages.
\fancypagestyle{plain}{%
    \fancyhf{}
    \renewcommand*\headrulewidth{0pt}
    \renewcommand*\footrulewidth{0pt}
    \fancyfoot[C]{\titlefont\thepage}
}

%% Fancy style for the main matter.
\setlength{\headheight}{14pt}
\fancypagestyle{mainmatter}{%
    \fancyhf{}
    %% Page numbers on the top left and top right.
    \fancyhead[LE,RO]{\titlefont\fontshape{n}\thepage}
    %% Chapter name on the left (even) page.
    \fancyhead[RE]{\titlefont\fontshape{n}\nouppercase{\leftmark}}
    %% Section name on the right (odd) page.
    \fancyhead[LO]{\titlefont\fontshape{n}\nouppercase{\rightmark}}
}

%% The mainmatter style is default for normal pages.
\pagestyle{mainmatter}

%% Print the current chapter and section at the top of the page in cyan.
\renewcommand*\chaptermark[1]{\markboth{\thechapter.\ \color{black}#1}{}}
\renewcommand*\sectionmark[1]{\markright{\thesection.\ \color{black}#1}}
\newcommand*\setheader[1]{\markboth{\color{title}#1}{\color{title}#1}}

%% Change the headrule command (from fancyhdr.sty) to draw the line below the
%% header in the title color.
\renewcommand*\headrule{%
    \if@fancyplain%
        \let\headrulewidth\plainheadrulewidth%
    \fi%
    {\color{title}\hrule\@height\headrulewidth\@width\headwidth}%
    \vskip-\headrulewidth%
}

%% Draw the line above a footnote in the title color as well.
\renewcommand*\footnoterule{%
    \vspace*{-3pt}%
    {\color{title}\hrule width 0.5\textwidth height 0.4pt}%
    \vspace*{2.6pt}%
}

%% A part title starts with a huge (96pt) bold black number, flushed to the
%% right, followed by the part name on the next line in the title color.
\titleformat{\part}[display]
    {\flushright}
    {\fontsize{96}{96}\selectfont\bfseries\thepart}
    {0pt}
    {\Huge\color{title}}
%% Separate the title from the text by two empty lines.
\titlespacing{\part}{0pt}{0pt}{2\baselineskip}
%% In the table of contents, the part name is preceded by an empty line, printed
%% in bold, and not followed by a line of dots.
\dottedcontents{part}[0em]{\vspace{\baselineskip}\titlefont\bfseries}{1.5em}{0pc}

%% Chapter titles have the same layout as parts.
\titleformat{\chapter}[display]
    {\flushright}
    {\firathin\Huge Chapter \firamedium\fontsize{96}{96}\selectfont\thechapter}
    {0pt}
    {\chaptitlestyle\Huge\color{title}}
\titlespacing{\chapter}{0pt}{0pt}{2\baselineskip}
%% In the table of contents, a chapter is similar to a part, except that it is
%% preceded by half an empty line.
\dottedcontents{chapter}[1.5em]{\vspace{0.5\baselineskip}\titlefont\bfseries}{1.5em}{0pc}

%% Section titles start with the number in bold, followed by the name printed
%% in the title color.
\titleformat{\section}
    {\Large\sectitlestyle}
    {\bfseries\firasemibold\thesection.\ }
    {0pt}
    {\color{title}}
%% Sections are preceded by an empty line.
\titlespacing{\section}{0pt}{\baselineskip}{0pt}
%% In the table of contents, section names are followed by a line of dots 8pt
%% apart.
\dottedcontents{section}[3.8em]{\titlefont}{2.3em}{8pt}

%% Subsection titles have the same layout as section titles, except in a smaller
%% font.
\titleformat{\subsection}
    {\large\sectitlestyle}
    {\bfseries\firasemibold\thesubsection.\ }
    {0pt}
    {\color{title}}
\titlespacing{\subsection}{0pt}{\baselineskip}{0pt}
\dottedcontents{subsection}[7em]{\titlefont}{3.2em}{8pt}

%% Subsubsections have the same font and color as sections and subsections, but
%% are not preceded by a number.
\titleformat{\subsubsection}
    {\headerstyle}
    {}
    {0pt}
    {\color{title}}
%% Subsubsections are preceded by an empty line and do not appear in the table
%% of contents.
\titlespacing{\subsubsection}{0pt}{\bigskipamount}{0pt}

%% Color the bullets of the itemize environment and make the symbol of the third
%% level a diamond instead of an asterisk.
\renewcommand*\labelitemi{\color{title}\textbullet}
\renewcommand*\labelitemii{\color{title}--}
\renewcommand*\labelitemiii{\color{title}$\diamond$}
\renewcommand*\labelitemiv{\color{title}\textperiodcentered}
