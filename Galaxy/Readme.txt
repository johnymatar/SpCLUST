Edit ~/galaxy/config/tool_conf.xml to add the modules (the edited tool_conf.xml present here)
Put the executables and sh files in a same directory (e.g. ~/SpCLUST)
Edit the XML files and fix the full path of the executable (e.g. /home/[your user directory]/SpCLUST/spclust-B62.sh $input $output)
Put the other XML files in ~/galaxy/tools/SpCLUST
Restart galaxy service if changes are made to tool_conf.xml
