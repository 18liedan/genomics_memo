# nig_supercomputer
Basic information about the NIG supercomputer and how to use it.
遺伝研スパコンの基本的な使い方

## General usage | 使い方
#### Steps on getting account, etc. | 利用開始までの流れはこちら
https://sc.ddbj.nig.ac.jp/guides/from_application_to_approval/
#### About nodes and logging in | ノードとログインについて
https://sc.ddbj.nig.ac.jp/guides/using_general_analysis_division/ga_usage/
In the general analysis division (free), you can use the "interactive node" and "compute node".
The interactive node can be used in the same way as you would use a regular Linux workstation.
You will use ssh to connect to it, and you can stay logged in for a maximum of 3 days.
一般区画（無料）ではインタラクティブノードと計算ノードが使える。
インタラクティブノードは、通常のLinuxワークスペースと同じようにリアルタイムで使える。
sshで繋ぐことができ、最大3日間ログインが継続される。
how to login into interactive node (e.g. my account name is "liedan"):
```
#open your terminal and type the following
ssh liedan@gw.ddbj.nig.ac.jp

#this should give you something that looks like this
	(base) liedan@gw: $

#then use ssh again to login to interactive node
ssh a001
	#you can also use a002, a003
```

For analyses that take longer than 3 days or need large amounts of CPU/RAM, you will use the compute node. I haven't yet used it myself, but in this node, you will schedule your jobs through SLURM.
計算ノードは、3日以上かかる解析や、CPUやメモリを大量に使う場合に使用する。まだ使ったことがないので詳しいことは解説できませんが、計算ノードではSLURMというジョブスケジューラーを使って、予約する。
About SLURM: https://sc.ddbj.nig.ac.jp/guides/software/JobScheduler/Slurm/

## General specs | 基本的なスペック
#### Storage space | ストレージについて：
Usually, the PI of the lab will be the representative of the group, and in total, one representative can get up to 30 TB for free, which will be allocated amongst account holders. (If you need more, you will need to pay). By default, each account will be allocated 1TB. If you would like to allocate more space per account, you need to apply.
通常、研究室のPIが責任者として複数のアカウント保持者を代表する。一責任者につき、30TBまで無料で使用でき、これをアカウント保持者に割り当てることができる。（30TB以上必要であれば、課金する。）特に何も申請しない場合、各アカウントにつき1TBが割り当てられる。それ以上が欲しい場合は、申請する。
Apply for more storage here: https://sc.ddbj.nig.ac.jp/application/resource_extension

Command for checking your storage (this example is for my "liedan" account):
`lfs quota -u liedan /home/liedan
(More details here: https://sc.ddbj.nig.ac.jp/guides/hardware/ga_lustre/)
#### CPU：
| Category                  | Value              | Explanation                                                     |
| ----------------------------- | ---------------------- | ------------------------------------------------------------------- |
| Physical Cores per Socket | 96                     | Each socket (CPU) has 96 physical cores.                            |
| Threads per Core          | 2                      | Each physical core supports 2 threads (hyperthreading enabled).     |
| Sockets                   | 2                      | There are 2 sockets (physical CPUs) installed in the system.        |
| Total Physical Cores      | 96×2=192               | Total physical cores across both sockets.                           |
| Max Logical Threads       | 192×2=384              | Total logical threads, accounting for hyperthreading-enabled cores. |
| CPU Frequency (Max)       | 3707.812 MHz (3.7 GHz) | Maximum supported clock speed per core.                             |

Command for checking CPU:
`lscpu`
```
#results from above command (2025 April 9th):
Architecture:             x86_64
  CPU op-mode(s):         32-bit, 64-bit
  Address sizes:          52 bits physical, 57 bits virtual
  Byte Order:             Little Endian

CPU(s):                   384
  On-line CPU(s) list:    0-383

Vendor ID:                AuthenticAMD
  Model name:             AMD EPYC 9654 96-Core Processor
    CPU family:           25
    Model:                17
    Thread(s) per core:   2
    Core(s) per socket:   96
    Socket(s):            2
    Stepping:             1
    Frequency boost:      enabled
    CPU(s) scaling MHz:   63%
    CPU max MHz:          3707.8120
    CPU min MHz:          1500.0000
    BogoMIPS:             4799.87
```
#### RAM | メモリ：
Command to check available RAM:
`free -h`
