- [Download AWS EC2 private key](#download-aws-ec2-private-key)
- [Access workshop EC2 server](#access-workshop-ec2-server)
  - [1. Connecting to EC2 workshop server using Putty](#1-connecting-to-ec2-workshop-server-using-putty)
  - [2. Connecting to Ubuntu desktop using RealVNC](#2-connecting-to-ubuntu-desktop-using-realvnc)
  - [(Optional) FileZilla client](#optional-filezilla-client)


## Download AWS EC2 private key

Please use [cumed_user_key.ppk](https://github.com/liu-xingliang/cumed_workshop/blob/main/00.prerequisites/keys/cumed_user_key.ppk) uploaded to workshop github repository (Rigt-click and save).

## Access workshop EC2 server

For X64 (64-bit x86) system (most user should use this one), please download: [putty](https://the.earth.li/~sgtatham/putty/latest/w64/putty.exe). For other architectures, please find corresponding binaries [here](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html).

### 1. Connecting to EC2 workshop server using Putty

1. Specify host name: `server_$(printf "%03d" N).cumed-workshop.com`, where `N` is a number between 1 and 30 (30 attendees), i.e., "server_001.cumed-workshop.com", "server_002.cumed-workshop.com", ..., "server_030.cumed-workshop.com"

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

5. Click "open" to connect (Please choose "connect once" in the subsequent window, so putty will not save current host keys, which will cause some problem in current workshop AWS settings: dynamic IP for the EC2 instance).

<p align="left">
<img src="./img/putty_connectonce.svg" width="400">
</p>

<p align="left">
<img src="./img/putty_terminal.svg" width="400">
</p>

### 2. Connecting to Ubuntu desktop using RealVNC

For x64 system, download [VNC-Viewer-7.8.0-Windows-64bit.exe](https://downloads.realvnc.com/download/file/viewer.files/VNC-Viewer-7.8.0-Windows-64bit.exe). Other architectures, please check the [download page](https://www.realvnc.com/en/connect/download/viewer/windows/). While connecting to the remote server using VNC through "localhost:5901", please **keep the putty ssh session alive**.

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

_PS: If attendees get VNC connection issues, please set up a fresh VNC session by running below commands on EC2 server:_

```bash
vncserver -kill :1
vncserver -localhost -geometry 1600x1200
```

### (Optional) FileZilla client

For x64 system, download [win64 installer](https://filezilla-project.org/download.php?platform=win64). For attendees use a lab PC **without** "Administrator" access:

1. Initiate the installation.

<p align="left">
<img src="./img/filezilla/1_click_run.svg" width="400">
</p>

2. Admin access is not required, please select "No".

<p align="left">
<img src="./img/filezilla/2_dont_need_admin.svg" width="400">
</p>

3. Only install for "me".

<p align="left">
<img src="./img/filezilla/3_agree_onlyforme.svg" width="800">
</p>

4. Specify a local folder to which current non-Admin user has access (e.g., Desktop).

<p align="left">
<img src="./img/filezilla/4_allcomp_choosefolder.svg" width="1000">
</p>

5. Setup a "New Site" (e.g., "cumed_workshop"):

* Use "SFTP" protocal.
* Specify host name.
* Logon type: "Key file".
* User: "ubuntu"
* Select .ppk key file we used for Putty ssh session.

<p align="left">
<img src="./img/filezilla/5_site.svg" width="800">
</p>

6. Please uncheck the box for NOT saving host key, as it might cause some problem in current workshop AWS settings (dynamic IP for the EC2 instance).

<p align="left">
<img src="./img/filezilla/6_uncheck_box.svg" width="400">
</p>

7. Direct to the target directory on the right panel and download (drag-and-drop is supported) selected files.

<p align="left">
<img src="./img/filezilla/7_direct_and_pull.svg" width="1000">
</p>









