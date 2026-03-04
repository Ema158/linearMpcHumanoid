#pragma once

enum class Contact {
    Both,
    Right,
    Left,
    None
};

class ContactState {
public:
    ContactState(Contact contact)
        : contactType_(contact) {}

    void set(Contact contact) {
        contactType_ = contact;
    }

    Contact get() const {
        return contactType_;
    }

    int getNumActiveContacts() const {
        switch (contactType_) {
            case Contact::Both:  return 2;
            case Contact::Right: return 1;
            case Contact::Left:  return 1;
            case Contact::None:  return 0;
        }
        return 0;
    }

    bool isRightActive() const {
        return contactType_ == Contact::Both ||
               contactType_ == Contact::Right;
    }

    bool isLeftActive() const {
        return contactType_ == Contact::Both ||
               contactType_ == Contact::Left;
    }

private:
    Contact contactType_;
};