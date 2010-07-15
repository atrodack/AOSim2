/*
 * Permission to use, copy, modify, and distribute this software and its
 * documentation under the terms of the GNU General Public License is hereby 
 * granted. No representations are made about the suitability of this software 
 * for any purpose. It is provided "as is" without express or implied warranty.
 * See the GNU General Public License for more details.
 *
 * Ian Martins (ianxm@jhu.edu)
 */

	int inCommentLine = 0;			/* in comment line */
	int comment_caller;
	void closeCommentBlock();
ID	[a-zA-Z][0-9a-zA-Z_]*
WS	[ \t]*

%s func funcDef paramList commentBlock
%%

<func>"%%"	{				/* end last function */
		  comment_caller = YY_START;
		  BEGIN(commentBlock);
		  fprintf(yyout,"}\n/**");
		  inCommentLine = 1;		/* start of comment */
		}
<commentBlock>"%%"		{
		  fprintf(yyout,"*/\n/**");	/* close old comment, start new one */
		  inCommentLine = 1;		/* start of comment */
		}
"%%"		{				/* first token */
		  comment_caller = YY_START;
		  BEGIN(commentBlock);
		  fprintf(yyout,"/**");
		  inCommentLine = 1;		/* start of comment */
		}
<commentBlock>"%"	{
		  if (inCommentLine)		/* just a symbol in a comment */
		    ECHO;
		  else 				/* part of current comment block */
		  {
		    fprintf(yyout,"*");
		    inCommentLine = 1;
		  }
		}
"%"		{				/* start comment block */
		  fprintf(yyout,"/*");
		  comment_caller = YY_START;
		  BEGIN(commentBlock);
		  inCommentLine = 1;		/* start of comment */
		}
\n		ECHO; inCommentLine = 0;
^{WS}function	{
		  closeCommentBlock(); 
		  fprintf(yyout,"void ");
		  BEGIN(funcDef);
		}
^{WS}function.*\[.*\].*=	{
		  closeCommentBlock(); 
		  fprintf(yyout,"rets ");
		  BEGIN(funcDef);
		}
^{WS}function.*=	{
		  closeCommentBlock(); 
		  fprintf(yyout,"ret ");
		  BEGIN(funcDef);
		}
<funcDef>\(	{
		  ECHO;
		  BEGIN(paramList);
		}
<paramList>\)	{
		  ECHO;
		  fprintf(yyout,"{\n");
		  BEGIN(func);
		}
<paramList>{ID}	{
		  fprintf(yyout,"type ");
		  ECHO;
		}
<commentBlock>[0-9a-zA-Z]+	{
		  closeCommentBlock();
		  ECHO;
		}
[0-9a-zA-Z]+	{
		  ECHO;
		}
<<EOF>>		{
		  fprintf(yyout,"}\n");
		  return;
		}

%%

void closeCommentBlock()			/* end comment block */
{
  if (!inCommentLine)
  {
    fprintf(yyout,"*/\n");
    BEGIN(comment_caller);
  }
}

main( int argc, char **argv )
{
  if ( argc > 1 )
    yyin = fopen( argv[1], "r" );
  else
    yyin = stdin;

  if ( argc > 2 )
    yyout = fopen( argv[2], "w" );
  else
    yyout = stdout;

  yylex();
}
