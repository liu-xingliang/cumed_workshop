## AWS EC2 instance and key pairs

This workshop will be conducted using [AWS EC2](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/concepts.html) instances. For each attendee, you will have your own EC2 instance with DNS name `server_$(printf "%03d" N).cumed-workshop.com`, where `N` is the rank on the registration list. E.g., if a attendee is ranked 1st on the list, the EC2 server assigned to him/her will be `server_001.cumed-workshop.com`. 

To access a AWS EC instance, attendee needs AWS [key pairs](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/create-key-pairs.html#having-ec2-create-your-key-pair), which has already been created. 

Attendee-specific EC2 server address and key pairs has already been sent to him/her through **registration email**.

## Access workshop EC2 server (VNC + Ubuntu Xfce desktop)

Workshop servers use Ubuntu 22.04 LTS with Xfce lightweight desktop env. After `ssh` to the server, attendees could use [VNC](https://en.wikipedia.org/wiki/Virtual_Network_Computing) to access the Xfce GUI desktop.

### Windows OS

For Windows users, please download putty/puttygen [binaries](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html):

For X64 (64-bit x86) system (most user should use this one):
>https://the.earth.li/~sgtatham/putty/latest/w64/putty.exe
>https://the.earth.li/~sgtatham/putty/latest/w64/puttygen.exe

For other architecture, please find corresponding binaries [here](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html).

#### 1. Coverting key pair (.pem) to Putty (.ppk) format

Double click puttygen binary downloaded and refer to section: **Convert your private key using PuTTYgen**.

Please follow the step-by-step instructions below:

1. Load .pem key.

<p align="left">
<img src="./img/putty_gen_load.png" width="400">
</p>

<p align="left">
<img src="./img/putty_load_2.png" width="400">
</p>

2. Save private key in .ppk format.

<p align="left">
<img src="./img/save_private_in_ppk.png" width="400">
</p>


#### 2. Connecting to EC2 workshop server using Putty

1. Specify host name.

<p align="left">
<img src="./img/putty_hostname.svg" width="400">
</p>

2. Left "Category" panel: "Connection" -> "SSH" -> "Auth" -> "Credentials"

<p align="left">
<img src="./img/putty_ppk.svg" width="400">
</p>

3. Enable port forwarding for VNC.

<p align="left">
<img src="./img/putty_portforwarding.svg" width="400">
</p>

4. Specify login user name.

<p align="left">
<img src="./img/putty_user.svg" width="400">
</p>

5. Click "open" to connect.

<p align="left">
<img src="./img/putty_terminal.svg" width="400">
</p>

#### 3. Connecting to Ubuntu desktop using RealVNC

Attendee could download a standalone binary of RealVNC client without installation (e.g., X64 system): https://downloads.realvnc.com/download/file/viewer.files/VNC-Viewer-7.8.0-Windows-64bit.exe. Other architectures, please check the [download page](https://www.realvnc.com/en/connect/download/viewer/windows/). While connecting to the remote server using VNC through "localhost:5901", please **keep the putty ssh session alive**.

<p align="left">
<img src="./img/windows_vnc.svg" width="400">
</p>

And click "continue" to proceed:

<p align="left">
<img src="./img/realvnc_continue.svg" width="400">
</p>

Password is "ubuntu":

<p align="left">
<img src="./img/realvnc_password.svg" width="400">
</p>

### Linux/MacOS

For Linux/MacOS users, please download the key pairs: `cumed_user_key.pem`and change permission:

```bash
curl -O https://raw.githubusercontent.com/liu-xingliang/cumed_workshop/main/00.prerequisites/cumed_user_key.pem
chmod 400 cumed_user_key.pem
```

After that, user could access their EC2 server with (e.g., for user 001), _please specify correct path of the key pair_:

```bash
ssh -i ~/Downloads/linux_EC2_key.pem -L 5901:localhost:5901 -o StrictHostKeyChecking=no ubuntu@aaa.cumed-workshop.com 
```

The `-L 5901:localhost:5901` part is to enable port forwarding to map local port 5901 to remote port 5901 upon which VNC server is listening. 

#### MacOS user VNC to EC2 server 

While keep the above ssh session alive. Attendees simply go "Finder" -> "Go" -> "Connect to Server" -> "localhost:5901" -> password: "ubuntu".

<p align="left">
<img src="./img/mac_vnc.svg" width="800">
</p>

#### Linux user VNC to EC2 server

Similar to Windows OS, attendees could also use standalone RealVNC client for Linux here: https://downloads.realvnc.com/download/file/viewer.files/VNC-Viewer-7.8.0-Linux-x64.

