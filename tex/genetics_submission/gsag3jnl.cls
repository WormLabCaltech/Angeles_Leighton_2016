% This class file is provided for authors submitting to the
% Genetics Society of America's G3: Genes | Genomes |
% Genetics Journal
%
% Version: 1.1, 2016/07/13
%
% Produced by Overleaf (www.overleaf.com)
%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{gsag3jnl}[13/07/2016, v1.6]
\AtEndOfClass{\RequirePackage[tracking=smallcaps]{microtype}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{extarticle}}
\ProcessOptions*
\LoadClass{extarticle}

% Fonts and language
\RequirePackage[utf8]{inputenc}
\RequirePackage[english]{babel}
\RequirePackage{amsmath,amsfonts,amssymb}
\RequirePackage{mathpazo}
\RequirePackage[scaled]{helvet}
\RequirePackage[T1]{fontenc}
\RequirePackage{url}
\RequirePackage[colorlinks=true, allcolors=blue]{hyperref}
\RequirePackage{pifont}

\newcommand{\watermarkfont}{\color{black50}\normalfont\fontsize{10}{12}\selectfont}
\newcommand{\headerfont}{\normalfont\sffamily\fontsize{8}{10}\selectfont}
\newcommand{\titlefont}{\color{color3}\normalfont\sffamily\bfseries\fontsize{21}{23}\selectfont}
\newcommand{\datesfont}{\normalfont\sffamily\itshape\fontsize{8}{10}\selectfont}
\renewcommand\UrlFont{\color{blue}\sffamily}
\newcommand{\absfont}{\normalfont\sffamily\linespread{1.2}\fontsize{10}{12}\selectfont}
\newcommand{\keywordstitlefont}{\color{color4}\normalfont\sffamily\linespread{1.2}\fontsize{8}{10}\selectfont}
\newcommand{\keywordsfont}{\normalfont\sffamily\linespread{1.2}\fontsize{10}{12}\selectfont}
\newcommand{\copyrightfont}{\normalfont\sffamily\fontsize{7}{9}\selectfont}
\newcommand{\footerfont}{\normalfont\sffamily\fontsize{8}{10}\selectfont}

\RequirePackage{datetime}

\newdateformat{monthyeardate}{%
  \monthname[\THEMONTH] \THEYEAR}

\RequirePackage{enumitem} 
\setlist[itemize]{noitemsep}


% Author and affiliation
\RequirePackage{authblk}
\setlength{\affilsep}{0.1em}
\renewcommand\Authfont{\color{color0}\normalfont\sffamily\bfseries\fontsize{9}{13}\selectfont}
\renewcommand\Affilfont{\color{color0}\normalfont\sffamily\fontsize{8}{10}\selectfont}
\makeatletter
\renewcommand\AB@affilsepx{, \protect\Affilfont}
\makeatother
\renewcommand\Authands{ and }

% Choose article type
\newcommand*{\articletype}[1]{%
  \RequirePackage{styles/#1}}
  
% Options for element switching
\RequirePackage{xifthen}
\newboolean{displaycopyright}

% Copyright statement
\RequirePackage{textcomp} % For copyright symbol styling
\makeatletter
\newcommand{\copyrightstatement}{Copyright \, \textcopyright\, \the\year{} \@correspondingauthor \hspace{0.2em} \textit{et al.} }
\makeatother

% Graphics, tables and other formatting
\RequirePackage{graphicx,xcolor}
\RequirePackage{colortbl}
\RequirePackage{booktabs}
\RequirePackage{tikz}
\RequirePackage{algorithm}
\RequirePackage[noend]{algpseudocode}
\RequirePackage{changepage}
\RequirePackage[left=36pt,%
                right=36pt,%
                top=48pt,%
                bottom=60pt,%
                headheight=40pt,%
                headsep=10pt,%
                letterpaper,twoside]{geometry}%
\RequirePackage[labelfont={bf,sf},%
                labelsep=space,%
                figurename=Figure,%
                singlelinecheck=off,%
                justification=RaggedRight]{caption}
\setlength{\columnsep}{24pt} % Distance between the two columns of text
\setlength{\parindent}{12pt} % Paragraph indent

% Set document color scheme
\definecolor{black50}{gray}{0.5} % 50% black for hrules
\definecolor{color0}{RGB}{0,0,0} % Base color
\definecolor{color1}{RGB}{246, 245, 209} % Table column heads background
% color2 defined by the specific article type style file
\definecolor{color3}{RGB}{0,55,123} % Article Titles
\definecolor{color4}{RGB}{134,0,55} % Subheadings: A heads and “Abstract”
\definecolor{tabledingbat}{RGB}{184,108,108} % Table dingbat icon

% Bibliography
\usepackage{natbib} \bibpunct{(}{)}{;}{author-year}{}{,}
\bibliographystyle{genetics}
\addto{\captionsenglish}{\renewcommand{\refname}{Literature Cited}}
\setlength{\bibsep}{0.0pt}

% Table style
\arrayrulecolor{color0} % Set the color of the horizontal rules
\renewcommand{\arraystretch}{1.5} % Increased line spacing
\captionsetup[table]{position=below} % Reduce spacing below caption
\colorlet{tableheadcolor}{color1} % Set table header colour
\newcommand{\header}{\rowcolor{tableheadcolor}} % create the command for the header
\captionsetup[table]{name={{\color{tabledingbat}\textrm{\ding{110}}}\hspace{0.8em}Table}} % Add the icon to the left of the table name

% Algorithm style
\makeatletter
\newcommand\fs@coloruled{\def\@fs@cfont{\bfseries}\let\@fs@capt\floatc@coloruled
  \def\@fs@pre{\kern2pt}%
  \def\@fs@post{\kern2pt{\color{color2}\hrule}\relax}%
  \def\@fs@mid{\kern2pt{\color{color2}\hrule}\kern2pt}%
  \let\@fs@iftopcapt\iftrue}
\makeatother

\floatstyle{coloruled}
\restylefloat{algorithm}

% Headers and footers
\RequirePackage{fancyhdr}  % custom headers/footers
\RequirePackage{lastpage}  % Number of pages in the document
\pagestyle{fancy}          % Enables the custom headers/footers

\fancypagestyle{firststyle}{
   \fancyhead[L]{}
   \fancyhead[R]{}
    \fancyfoot[L]{\includegraphics[width=96pt]{G3_Horizontal_Color_High_Resolution.jpg}}
   \fancyfoot[R]{{\footerfont Volume X\hspace{1em}|\hspace{1em}\monthyeardate\today\hspace{1em}|\hspace{1em}\thepage\ }}
}

% Headers
\fancyhead[LO,RE,LE,RO]{} % Blank headers for non-title pages.
% Footers
\lfoot{}
\cfoot{}
\rfoot{}
\makeatletter
\fancyfoot[LE]{{\footerfont \thepage\hspace{1em}|\hspace{1em}\@correspondingauthor \hspace{0.2em} \textit{et al.} }}
\fancyfoot[RO]{{\footerfont Volume X\hspace{1em}\monthyeardate\today\hspace{1em}|\hspace{1em}\@runningtitle\hspace{1em}|\hspace{1em}\thepage\ }}
\makeatother
\fancyfoot[LO,RE]{\includegraphics[width=96pt]{G3_Horizontal_Color_High_Resolution.jpg}}
\renewcommand{\headrulewidth}{0pt}% % No header rule
\renewcommand{\footrulewidth}{0pt}% % No footer rule

% Article type header box
\newcommand{\logomark}{
\begin{textblock}{200}(36,36)
\begin{flushleft}
\includegraphics[width=126pt]{G3_Horizontal_Color_High_Resolution.jpg}
\end{flushleft}
\end{textblock}
}

% Article type header box
\newcommand{\articletypemark}{
\begin{textblock}{200}(378,36)
\begin{flushright}
{\raggedleft \headerfont\color{black} \textls[96]{\articletypename}}
\end{flushright}
\end{textblock}
}

% Colored margin box
\newcommand{\marginmark}{
\begin{textblock}{72}(506,0) 
\noindent
\begin{tikzpicture}[x=1,y=1]
\fill[color=color2] (0,0) rectangle (72,36);
\end{tikzpicture}
\end{textblock}
}

%% 13 July 2016: make sure sections are unnumbered
\setcounter{secnumdepth}{0}
\setcounter{tocdepth}{0}

% Section/subsection/paragraph set-up
\RequirePackage[explicit]{titlesec}
\renewcommand{\thesubsection}{\Alph{subsection}}

\titleformat{\section}
  {\large\sffamily\bfseries\color{color4}}
  {\thesection.}
  {0.5em}
  {\MakeUppercase{#1}}
  []
\titleformat{name=\section,numberless}
  {\large\sffamily\bfseries\color{color4}}
  {}
  {0em}
  {\MakeUppercase{#1}}
  []
\titleformat{\subsection}
  {\sffamily\bfseries}
  {\thesubsection.}
  {0.5em}
  {#1}
  []
\titleformat{\subsubsection}[runin]
  {\sffamily\bfseries\itshape}
  {\thesubsubsection.}
  {0.5em}
  {#1}
  []
\titleformat{\paragraph}[runin]
  {\sffamily\small\bfseries}
  {}
  {0em}
  {#1}
\titlespacing*{\section}{0pc}{3ex \@plus4pt \@minus3pt}{5pt}
\titlespacing*{\subsection}{0pc}{2.5ex \@plus3pt \@minus2pt}{2pt}
\titlespacing*{\subsubsection}{0pc}{2ex \@plus2.5pt \@minus1.5pt}{2pt}
\titlespacing*{\paragraph}{0pc}{1.5ex \@plus2pt \@minus1pt}{12pt}

% Article meta data additional fields
\newcommand{\correspondingauthor}[1]{\def\@correspondingauthor{#1}}
\newcommand{\keywords}[1]{\def\@keywords{#1}}
\newcommand{\runningtitle}[1]{\def\@runningtitle{#1}}

\newcommand{\keywordname}{KEYWORDS} % Defines the keywords heading name

\def\xabstract{abstract}
\long\def\abstract#1\end#2{\def\two{#2}\ifx\two\xabstract
\long\gdef\theabstract{\ignorespaces#1}
\def\go{\end{abstract}}\else
\typeout{^^J^^J PLEASE DO NOT USE ANY \string\begin\space \string\end^^J
COMMANDS WITHIN ABSTRACT^^J^^J}#1\end{#2}
\gdef\theabstract{\vskip12pt BADLY FORMED ABSTRACT: PLEASE DO
NOT USE {\tt\string\begin...\string\end} COMMANDS WITHIN
THE ABSTRACT\vskip12pt}\let\go\relax\fi
\go}

% Set environment for abstract content and styling
\newcommand{\abscontent}{
\begin{minipage}[t]{453pt}
{\color{black50}\rule{\linewidth}{0.5pt}}
\vskip9pt%
{\absfont\bfseries\color{color4} ABSTRACT} \absfont \theabstract%
\vskip1pt%
{\color{black50}\rule{\linewidth}{0.5pt}}
\end{minipage}\hspace{15pt}
\begin{minipage}[t]{72pt}
\vskip18.3pt%
{\keywordstitlefont\bfseries\keywordname} \\ \keywordsfont \@keywords%
\end{minipage}
}%

% Custom title page
\renewcommand{\@maketitle}{\bgroup\setlength{\parindent}{0pt}
\begin{adjustwidth}{0pt}{24pt}
\begin{flushleft}
{%
\vskip96pt%
{\raggedright \titlefont \@title\par}%
\vskip11pt
{\raggedright \@author\par}
\vskip20pt%
}%
\end{flushleft}
\end{adjustwidth}
\egroup
{%
\abscontent
}%
\vskip25pt%
}%

% Other packages
\RequirePackage{enumitem} % For reducing bullet list item separation

\RequirePackage{tabularx}

\RequirePackage[flushmargin,ragged]{footmisc}

\newlength{\myFootnoteWidth}%
\newlength{\myFootnoteLabel}%
\setlength{\myFootnoteLabel}{0.5em}
\renewcommand{\@makefntext}[1]{%
  \setlength{\myFootnoteWidth}{\columnwidth}%
  \addtolength{\myFootnoteWidth}{-\myFootnoteLabel}%
  \noindent\makebox[\myFootnoteLabel][r]{\@makefnmark\ }%
  \parbox[t]{\myFootnoteWidth}{#1}
  }

% Commands to set the footnote for the front page.
\renewcommand{\footnoterule}{%
  \kern -3pt
  {\color{black50} \hrule width \columnwidth height 0.5pt}
  \kern 2.5pt
}

\newcommand\blfootnote[1]{%
  \begingroup
  \setlength{\myFootnoteLabel}{0em}
  \renewcommand\thefootnote{}\footnote{#1}%
  \addtocounter{footnote}{-1}%
  \addtocounter{Hfootnote}{-1}%
  \endgroup
}

\newcommand\correspondingauthoraffiliation[1]{%
  \begingroup
  \setlength{\myFootnoteLabel}{0em}
\renewcommand\thefootnote{}\footnote{\footnotelayout\copyrightfont\textsuperscript{1}#1}%
  \addtocounter{footnote}{-1}%
  \addtocounter{Hfootnote}{-1}%
  \endgroup
}

\newcommand{\firstpagefootnote}{\blfootnote{\copyrightfont \ifthenelse{\boolean{displaycopyright}}{\copyrightstatement}{}\\%
Manuscript compiled: \today\%
}}

\newenvironment{tableminipage}[1]{\begin{minipage}{#1}\renewcommand\footnoterule{ \kern -1ex}%
\setlength{\myFootnoteLabel}{0.5em}%
}{\end{minipage}}

% Package to enable absolute positioning of figures
\RequirePackage[absolute,overlay]{textpos}
\setlength{\TPHorizModule}{1pt}
\setlength{\TPVertModule}{1pt}
\textblockorigin{0pt}{0pt} % start everything near the top-left corner