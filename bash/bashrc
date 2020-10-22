#
[[ $- != *i* ]] && return

colors() {
	local fgc bgc vals seq0

	printf "Color escapes are %s\n" '\e[${value};...;${value}m'
	printf "Values 30..37 are \e[33mforeground colors\e[m\n"
	printf "Values 40..47 are \e[43mbackground colors\e[m\n"
	printf "Value  1 gives a  \e[1mbold-faced look\e[m\n\n"

	# foreground colors
	for fgc in {30..37}; do
		# background colors
		for bgc in {40..47}; do
			fgc=${fgc#37} # white
			bgc=${bgc#40} # black

			vals="${fgc:+$fgc;}${bgc}"
			vals=${vals%%;}

			seq0="${vals:+\e[${vals}m}"
			printf "  %-9s" "${seq0:-(default)}"
			printf " ${seq0}TEXT\e[m"
			printf " \e[${vals:+${vals+$vals;}}1mBOLD\e[m"
		done
		echo; echo
	done
}

[ -r /usr/share/bash-completion/bash_completion ] && . /usr/share/bash-completion/bash_completion

# Change the window title of X terminals
case ${TERM} in
	xterm*|rxvt*|Eterm*|aterm|kterm|gnome*|interix|konsole*)
		PROMPT_COMMAND='echo -ne "\033]0;${USER}@${HOSTNAME%%.*}:${PWD/#$HOME/\~}\007"'
		;;
	screen*)
		PROMPT_COMMAND='echo -ne "\033_${USER}@${HOSTNAME%%.*}:${PWD/#$HOME/\~}\033\\"'
		;;
esac

use_color=true

# Set colorful PS1 only on colorful terminals.
# dircolors --print-database uses its own built-in database
# instead of using /etc/DIR_COLORS.  Try to use the external file
# first to take advantage of user additions.  Use internal bash
# globbing instead of external grep binary.
safe_term=${TERM//[^[:alnum:]]/?}   # sanitize TERM
match_lhs=""
[[ -f ~/.dir_colors   ]] && match_lhs="${match_lhs}$(<~/.dir_colors)"
[[ -f /etc/DIR_COLORS ]] && match_lhs="${match_lhs}$(</etc/DIR_COLORS)"
[[ -z ${match_lhs}    ]] \
	&& type -P dircolors >/dev/null \
	&& match_lhs=$(dircolors --print-database)
[[ $'\n'${match_lhs} == *$'\n'"TERM "${safe_term}* ]] && use_color=true

if ${use_color} ; then
	# Enable colors for ls, etc.  Prefer ~/.dir_colors #64489
	if type -P dircolors >/dev/null ; then
		if [[ -f ~/.dir_colors ]] ; then
			eval $(dircolors -b ~/.dir_colors)
		elif [[ -f /etc/DIR_COLORS ]] ; then
			eval $(dircolors -b /etc/DIR_COLORS)
		fi
	fi

	if [[ ${EUID} == 0 ]] ; then
		PS1='\[\033[01;31m\][\h\[\033[01;36m\] \W\[\033[01;31m\]]\$\[\033[00m\] '
	else
		PS1='\[\033[01;32m\][\u@\h\[\033[01;37m\] \W\[\033[01;32m\]]\$\[\033[00m\] '
	fi

	alias ls='ls --color=auto'
	alias grep='grep --colour=auto'
	alias egrep='egrep --colour=auto'
	alias fgrep='fgrep --colour=auto'
else
	if [[ ${EUID} == 0 ]] ; then
		# show root@ when we don't have colors
		PS1='\u@\h \W \$ '
	else
		PS1='\u@\h \w \$ '
	fi
fi

unset use_color safe_term match_lhs sh

alias cp="cp -i"                          # confirm before overwriting something
alias df='df -h'                          # human-readable sizes
alias free='free -m'                      # show sizes in MB
alias np='nano -w PKGBUILD'
alias more=less

xhost +local:root > /dev/null 2>&1

complete -cf sudo

# Bash won't get SIGWINCH if another process is in the foreground.
# Enable checkwinsize so that bash will check the terminal size when
# it regains control.  #65623
# http://cnswww.cns.cwru.edu/~chet/bash/FAQ (E11)
shopt -s checkwinsize

shopt -s expand_aliases

# export QT_SELECT=4

# Enable history appending instead of overwriting.  #139609
shopt -s histappend

#
# # ex - archive extractor
# # usage: ex <file>
ex ()
{
  if [ -f $1 ] ; then
    case $1 in
      *.tar.bz2)   tar xvjf $1   ;;
      *.tar.gz)    tar xvzf $1   ;;
      *.bz2)       bunzip2 $1   ;;
      *.rar)       unrar x $1     ;;
      *.gz)        gunzip $1    ;;
      *.tar)       tar xvf $1    ;;
      *.tbz2)      tar xvjf $1   ;;
      *.tgz)       tar xvzf $1   ;;
      *.zip)       unzip $1     ;;
      *.Z)         uncompress $1;;
      *.7z)        7z x $1      ;;
      *)           echo "'$1' cannot be extracted via ex()" ;;
    esac
  else
    echo "'$1' is not a valid file"
  fi
}


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/angel/.anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/angel/.anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/angel/.anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/angel/.anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<



# -- ADDED -----------------------------------------------------------

# - To unpack downloaded files
funpack_all ()
{
  for i in `ls | grep fits | grep hmi`
  do 
    mv $i $i.fz
    funpack -E 1 -v $i.fz
    rm $i.fz
  done
}


# -- To propertly compile PDF files from .tex sources
pdf ()
{
  filename=$(basename $1)
  filename="${filename%.*}"
  if [ -f $1 ];
  then
    for i in {1..3};
    do
      pdflatex $filename.tex
      bibtex $filename.aux 
    done
  fi
  rm *MSc.aux *.bbl *.blg *.lof *MSc.log *.out *.synctex.gz *.toc 
  mv $filename.pdf out.$filename.pdf
  evince out.$filename.pdf 
}

#alias ls='ls -lh --color --group-directories-first'
alias la='ls -lh --color --group-directories-first'
alias l='ls -lh --color --group-directories-first'
alias ..='cd ..'
alias ...='cd ../..'
alias ....='cd ../../..'

alias master='cd $HOME/Documents/MASTER'
alias masterc='cd $HOME/Documents/MASTER/Codes'
alias mastercc='cd $HOME/Documents/MASTER/Codes/Charlie'
alias mastercj='cd $HOME/Documents/MASTER/Codes/JOBS'
alias mastercb='cd $HOME/Documents/MASTER/Codes/bash'
alias mastercp='cd $HOME/Documents/MASTER/Codes/Python'
alias flare='cd $HOME/IDLWorkspace/Flare'
alias prename='perl-rename'
#alias ts='texstudio $1; rm *.toc'

# -- Open TexStudio
tsd ()
{
  filename=$(basename $1)
  filename="${filename%.*}"
  texstudio $filename.tex
  rm *MSc.aux *.bbl *.blg *.lof *MSc.log *.out *.synctex.gz *.toc 
  mv $filename.pdf out.$filename.pdf
}

alias pfm=pcmanfm

for i in {1..69}
do
  alias 'f'$i='cd $HOME/IDLWorkspace/Flare/F'$i
  alias 'f'$i'd'='cd $HOME/IDLWorkspace/Flare/F'$i'/HMIDoppler'
  alias 'f'$i'dw'='cd $HOME/IDLWorkspace/Flare/F'$i'/HMIDoppler/Work'
  alias 'f'$i'i'='cd $HOME/IDLWorkspace/Flare/F'$i'/HMIIntensity'
  alias 'f'$i'iw'='cd $HOME/IDLWorkspace/Flare/F'$i'/HMIIntensity/Work'
  alias 'f'$i'm'='cd $HOME/IDLWorkspace/Flare/F'$i'/HMIMagneto'
  alias 'f'$i'mw'='cd $HOME/IDLWorkspace/Flare/F'$i'/HMIMagneto/Work'
  alias 'pf'$i'd'='pfm $HOME/IDLWorkspace/Flare/F'$i'/HMIDoppler'
  alias 'pf'$i'dw'='pfm $HOME/IDLWorkspace/Flare/F'$i'/HMIDoppler/Work'
  alias 'pf'$i'iw'='pfm $HOME/IDLWorkspace/Flare/F'$i'/HMIIntensity/Work'
  alias 'pf'$i'mw'='pfm $HOME/IDLWorkspace/Flare/F'$i'/HMIMagneto/Work'
done


# added by Angel for Charlie codes
export PATH="$HOME/Documents/MASTER/Codes/Charlie:$PATH"
export PATH="$HOME/Documents/MASTER/Codes/Charlie/shellprogs:$PATH"
export PATH="$HOME/Documents/MASTER/Codes/JOBS:$PATH"
export PATH="$HOME/Documents/MASTER/Codes/GitHub/JOBS:$PATH"
export PATH="$HOME/Documents/MASTER/Codes/Python:$PATH"
export PATH="$HOME/Documents/MASTER/Codes/Python/TD:$PATH"
#set PATH=($HOME/bin $path)

export PYTHONPATH=$HOME/Downloads/PLUTO/Tools/pyPLUTO/pyPLUTO
export PYTHONPATH=$PYTHONPATH:$HOME/ME-I_Workshop/pyMilne
export PYTHONPATH=$PYTHONPATH:$HOME/Documents/MASTER/Codes/GitHub/config_files/Python
export PYTHONPATH=$PYTHONPATH:$HOME/Documents/MASTER/Codes/Python
#export PYTHONPATH

export GITHUB_DIR=$HOME"/Documents/MASTER/Codes/GitHub" 
export GITPY_DIR=$HOME"/Documents/MASTER/Codes/GitHub/config_files/Python"

# -- Overall run of cfitsio programs (like funpack/fpack)
export LD_LIBRARY_PATH=$HOME/Documents/MASTER/Codes/cfitsio-3.49:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$HOME/Documents/MASTER/Codes/cfitsio-3.49/lib/pkgconfig:$LD_LIBRARY_PATH
export PATH="$HOME/Documents/MASTER/Codes/cfitsio-3.49:$PATH"
#export PATH="$HOME/Documents/MASTER/Codes/cfitsio-3.49/lib/pkgconfig:$PATH"
export PATH="$HOME/Documents/MASTER/Codes/bash:$PATH"

# -- To no repeat commands on bash history
export HISTCONTROL=ignoreboth:erasedups

# -- Add to proper latex compilation (put .sty files here)
export TEXMFHOME=$HOME/Documents/MASTER/texmf

HISTFILESIZE=20000
HISTSIZE=20000
export JAVA_HOME=/usr/lib/jvm/java-11-openjdk/
export PLUTO_DIR="$HOME/Downloads/PLUTO"
PROMPT_COMMAND=${PROMPT_COMMAND:+$PROMPT_COMMAND; }'printf "\033]0;%s\007" "${PWD/#$HOME/\~}"'
#PROMPT_COMMAND=${PROMPT_COMMAND:+$PROMPT_COMMAND; }'printf "\033]0;%s@%s:%s\007" "${USER}" "${HOSTNAME%%.*}" "${PWD/#$HOME/\~}"'
