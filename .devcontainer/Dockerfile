FROM ubuntu:latest

RUN apt update -y
RUN apt upgrade -y

RUN apt install git -y
RUN apt install g++ -y
RUN apt install make -y
RUN apt install zsh -y
RUN apt install sudo -y
RUN apt install clang-format -y



RUN useradd --create-home --password vscode --shell /bin/zsh vscode
COPY ./.zshrc /home/vscode


