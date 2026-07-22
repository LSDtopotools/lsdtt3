# Notes on git from the command line for new users

These are notes for people unfamiliar with git and GitHub who need to clone the
repository — for beta testing, or to build from source.

GitHub no longer accepts username/password authentication over HTTPS, so you
need an SSH key. These notes go through it step by step.

> **This is written for a Linux environment.** On Windows use WSL2, and on
> macOS use an Ubuntu VM; see
> [`building_from_source.md`](building_from_source.md) for both.

---

## Clone using an SSH key

### 1. Generate an SSH key pair

Skip this if you already have one — check with `ls ~/.ssh/*.pub`.

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
# Press Enter to accept the default file location
# Add a passphrase for security
```

Replace `your_email@example.com` with your actual email address. The key pair
is saved to the default location: `~/.ssh/id_ed25519` (the **private** key,
never share it) and `~/.ssh/id_ed25519.pub` (the **public** key).

> Older guides use `-t rsa`, which produces `~/.ssh/id_rsa`. Ed25519 is the
> current recommendation — shorter, faster and at least as secure. If you have
> an existing `id_rsa` key it still works; just substitute the filename
> throughout.

### 2. Add the key to the SSH agent

The agent holds the decrypted key so you are not asked for the passphrase on
every operation. Start it first if it is not already running:

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
```

The agent does not survive a reboot. If you get asked for your passphrase
repeatedly in later sessions, run these two lines again.

### 3. Copy your public key to GitHub

Print the **public** key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Copy the entire line, including the `ssh-ed25519` prefix and the trailing
comment. Then:

1. Go to your GitHub profile settings.
2. Click **SSH and GPG keys**.
3. Click **New SSH key**.
4. Paste the key into the **Key** field.
5. Give it a title (e.g. "Work laptop" — something that tells you which machine
   it is, so you can revoke the right one later).
6. Click **Add SSH key**.

### 4. Check that it works

```bash
ssh -T git@github.com
```

Expect: `Hi <username>! You've successfully authenticated, but GitHub does not
provide shell access.` That message is a success, despite how it reads.

The first connection asks you to confirm GitHub's host fingerprint — answer
`yes`.

### 5. Clone the repository

```bash
git clone git@github.com:lsdtopotools/lsdtt3.git
```

More generally, replace the owner and repository name as needed:

```bash
git clone git@github.com:your_username/repository_name.git
```

Note the `git@github.com:` form — an SSH URL. The `https://github.com/...`
form is what triggers the password prompt that no longer works.

---

## Fixing a repository you already cloned over HTTPS

If you cloned before setting up SSH, the remote is still an HTTPS URL and will
keep asking for a password. Point it at the SSH URL instead:

```bash
git remote -v                                              # check what it is now
git remote set-url origin git@github.com:lsdtopotools/lsdtt3.git
git remote -v                                              # confirm the change
```

The subcommand is `git remote set-url`, not `git set-url`.

---

## See also

- [`building_from_source.md`](building_from_source.md) — what to do once you have cloned
- [`developer_onboarding.md`](developer_onboarding.md) — getting started as a contributor
