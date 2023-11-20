## AWS EC2 instance and key pairs

This workshop will be conducted using [AWS EC2](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/concepts.html) instances. For each attendee, you will have your own EC2 instance with DNS name `server_$(printf "%03d" N).cumed-workshop.com`, where `N` is the rank on the registration list. E.g., if a attendee is ranked 1st on the list, the EC2 server assigned to him/her will be `server_001.cumed-workshop.com`. 

To access a AWS EC instance, attendee needs AWS [key pairs](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/create-key-pairs.html#having-ec2-create-your-key-pair), which has already been created. 

Attendee-specific EC2 server address and key pairs has already been sent to him/her through **registration email**.

## Access workshop EC2 server

### Linux/MacOS

For Linux/MacOS users, please download the key pairs: `cumed_user_key.pem_1`. Rename it and change permission:

```bash
mv cumed_user_key.pem_1 cumed_user_key.pem
chmod 400 cumed_user_key.pem
```

After that, user could access their EC2 server with (e.g., for user 001), _please specify correct path of the key pair_:

```bash
ssh -i /path/to/cumed_user_key.pem ubuntu@server_001.cumed-workshop.com 
```

For details, please refer [here](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/connect-linux-inst-ssh.html)

### Windows

For Windows users, please download putty/puttygen binaries:

For X64 (64-bit x86) system (most user should use this one):
>https://the.earth.li/~sgtatham/putty/latest/w64/putty.exe
>https://the.earth.li/~sgtatham/putty/latest/w64/puttygen.exe

For X86 system (32-bit x86) system:
>https://the.earth.li/~sgtatham/putty/latest/w32/putty.exe
>https://the.earth.li/~sgtatham/putty/latest/w32/puttygen.exe

For other architecture, please find corresponding binaries [here](https://www.chiark.greenend.org.uk/~sgtatham/putty/latest.html).

#### 1. Coverting key pair (.pem) to Putty (.ppk) format

Double click puttygen binary downloaded and refer to section: **Convert your private key using PuTTYgen** of [AWS guide](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/putty.html#putty-prereqs).

#### 2. Connecting to EC2 workshop server using Putty

Use **convered** key pair and refer to [**Connect to your Linux instance**](https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/putty.html#putty-ssh).

