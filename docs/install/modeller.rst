*******************************
Modeller --- Homology Modelling
*******************************

.. highlight:: bash


Installation procedure
======================

For the impatients
------------------

Open a shell and type (v9.17+)::

    version=9.17
    mkdir $HOME/Downloads/modeller; cd $_
    wget -O - http://salilab.org/modeller/${version}/modeller-${version}.tar.gz | tar xfz - -C . --strip-components=1
    ./Install <<EOF
    
    ~/bin/modeller${version}
    <password>
    
    
    EOF
    cd ..; rm -rf modeller

Sourcing in your .bashrc file::

    # added by $(whoami) on $(date) to source MODELLER
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:~/bin/modeller${version}/lib/x86_64-intel8"
    export PYTHONPATH="$PYTHONPATH:~/bin/modeller${version}/modlib:~/bin/modeller${version}/lib/x86_64-intel8/python2.5"

Links (< v9.16)::

    # shortcuts
    ln -sf ~/bin/modeller${version}/bin/mod${version} ~/bin/modeller
    ln -sf ~/bin/modeller${version}/bin/mod${version} ~/bin/mod${version}
    # .bashrc
    echo -e "\n# added by $(whoami) on $(date) to source MODELLER libraries" >> ~/.bashrc
    echo export LD_LIBRARY_PATH=\"$(find ~/bin/modeller${version}/lib/ -maxdepth 2 -name _modeller.so -printf '%h'):\$LD_LIBRARY_PATH\" >> ~/.bashrc
    # for Python2.7
    ln -sf ~/bin/modeller${version}/modlib/modeller              ~/.local/lib/python2.7/site-packages/modeller
    ln -sf ~/bin/modeller${version}/lib/*/python2.5/_modeller.so ~/.local/lib/python2.7/site-packages/_modeller.so
    # for Anaconda
    ln -sf ~/bin/modeller${version}/modlib/modeller              ~/.anaconda/lib/python2.7/site-packages/modeller
    ln -sf ~/bin/modeller${version}/lib/*/python2.5/_modeller.so ~/.anaconda/lib/python2.7/site-packages/_modeller.so



