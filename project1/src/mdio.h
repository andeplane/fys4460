#pragma once
class System;
class Settings;

class MDIO
{
public:
    System *system;
    Settings *settings;
    MDIO();
    void setup(System *system_);
    void save_state_to_file_binary();
};
