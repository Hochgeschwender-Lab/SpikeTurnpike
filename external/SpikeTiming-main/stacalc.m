% stacalc               compute multi-sta from tagged spike trains and an analog vector
%
% call                  [ f, fsem ] = stacalc( t, x, win, tag, utags )
%
% gets                  t, x, win   see sta.m
%                       tag         same dimensions as t; integers
%                       utags       the tags to be considered
%
% calls                 sta, stwiener, wsta
%
% Note: t can also be a tagged sparse matrix; then the tag argument is ignored 
%
% see also              isisort, stahat

% since the simple STA filter ignores history/other cells, another option is the
% cSTA filter: this is akin to the "response-conditional ensembles" of Bialek
% (de Ruyter van Steveninck and Bialek, 1988; Rieke et al., 1997; 
% see also Sharpee et al., 2006).

% 27-oct-14 ES

% last update
% September 2021

function [ f, fsem ] = stacalc( t, x, win, tag, utags, method )

% arguments
nargs                               = nargin;
if nargs < 3 || isempty( t ) || isempty( x ) || isempty( win )
    f                               = [];
    fsem                            = [];
    return
end
if length( win ) ~= 2 || ~isequal( win, sort( win ) )
    error( 'win must be two element vector' )
end 
nwin                                = diff( win ) + 1;
if issparse( t ) 
    tvec                            = find( t );
    tag                             = full( t( tvec ) );
    t                               = tvec;
end
if nargs < 4 || isempty( tag )
    tag                             = ones( size( t ) );
end
if ~isequal( size( t ), size( tag ) )
    error( 'input size mismatch' )
end
if nargs < 5 || isempty( utags )
    utags                           = unique( tag );
end
ntags                               = length( utags );
if nargs < 6 || isempty( method ) || ~isa( method, 'char' )
    method                          = 'sta';
end
method                              = lower( method );

% initialization
f                                   = NaN * ones( nwin, ntags );
fsem                                = f;

% sta
for i                               = utags( : ).'
    idx                             = tag == i;
    switch method
        case 'sta'
            [ f( :, i ), ss, nn ]   = sta( t( idx ), x, win );
            fsem( :, i )            = ss / sqrt( nn );
        case 'wsta'
            [ f( :, i ), fsem( :, i ) ] = wsta( t( idx ), x, win );
        case 'wiener'
            [ f( :, i ), fsem( :, i ) ] = stwiener( t( idx ), x, win );
    end
end

return

% EOF